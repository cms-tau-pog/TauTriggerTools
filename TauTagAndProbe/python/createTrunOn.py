#!/usr/bin/env python

import argparse
from array import array
import math
import numpy as np
import os
import re
import sys
import ROOT

parser = argparse.ArgumentParser(description='Create turn on curves.')
parser.add_argument('--input-data', required=True, type=str, help="skimmed data input")
parser.add_argument('--input-dy-mc', required=True, type=str, help="skimmed DY MC input")
parser.add_argument('--output', required=True, type=str, help="output file prefix")
parser.add_argument('--channels', required=False, type=str, default='etau,mutau,ditau', help="channels to process")
parser.add_argument('--decay-modes', required=False, type=str, default='all,0,1,10,11', help="decay modes to process")
parser.add_argument('--working-points', required=False, type=str,
                    default='VVVLoose,VVLoose,VLoose,Loose,Medium,Tight,VTight,VVTight',
                    help="working points to process")
parser.add_argument('--branchname-weight-data', required=True, type=str, help="branchname for event weights for data input")
parser.add_argument('--branchname-weight-dy-mc', required=True, type=str, help="branchname for event weights for DY MC input")

args = parser.parse_args()

if not(args.branchname_weight_data == "weight" or args.branchname_weight_data == "final_weight"):
    raise ValueError("Invalid configuration parameter branchname-weight-data = '%s' !!" % args.branchname_weight_data)
if not(args.branchname_weight_dy_mc == "weight" or args.branchname_weight_dy_mc == "final_weight"):
    raise ValueError("Invalid configuration parameter branchname-weight-dy-mc = '%s' !!" % args.branchname_weight_dy_mc)

path_prefix = '' if 'TauTriggerTools' in os.getcwd() else 'TauTriggerTools/'
sys.path.insert(0, path_prefix + 'Common/python')
from AnalysisTypes import *
from AnalysisTools import *
import RootPlotting
ROOT.ROOT.EnableImplicitMT(4)
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2() # CV: This does not seem to work... all bin-error are the exact square-root of the bin-contents !!
RootPlotting.ApplyDefaultGlobalStyle()

def CreateBins(max_pt, for_fitting):
    bins = None
    epsilon = 1.e-3 # CV: increase x-axis range a little bit, to make sure that upper limit (120 GeV) gets added to list of bin-edges
    if for_fitting:
        bins = np.arange(20, 40, step=1) 
        bins = np.append(bins, np.arange(40, 60, step=2))
        bins = np.append(bins, np.arange(60, 80, step=5))
        bins = np.append(bins, np.arange(80, 120, step=10))
        bins = np.append(bins, np.arange(120, 160, step=20))
        bins = np.append(bins, np.arange(160, 200+epsilon, step=40))
    else:
        bins = np.arange(20, 40, step=1)
        bins = np.append(bins, np.arange(40, 60, step=2))
        bins = np.append(bins, np.arange(60, 80, step=5))
        bins = np.append(bins, np.arange(80, 120, step=10))
        bins = np.append(bins, np.arange(120, 160, step=20))
        bins = np.append(bins, np.arange(160, 200+epsilon, step=40))
    use_logx = max_pt > 200+epsilon
    return bins, use_logx

class TurnOnData:
    def __init__(self):
        self.hist_total = None
        self.hist_passed = None
        self.eff = None

##def dumpHistogram(histogram):
##    for idxBin in range(histogram.GetNbinsX()):
##        print(" bin #%i: bin-content = %1.2f +/- %1.2f" % (idxBin + 1, histogram.GetBinContent(idxBin + 1), histogram.GetBinError(idxBin + 1)))

def CreateHistograms(input_file, branchname_weight,
                     channels, decay_modes, discr_name, 
                     working_points, hist_models, label, var):
    ##print("<CreateHistograms>:")
    ##print(" input_file = '%s'" % input_file)
    ##print(" branchname_weight = '%s'" % branchname_weight)
    df = ROOT.RDataFrame('events', input_file)
    turnOn_data = {}
    dm_labels = {}

    for dm in decay_modes:
        if dm == 'all':
            dm_labels[dm] = ''
            df_dm = df
        else:
            dm_labels[dm] = '_dm{}'.format(dm)
            df_dm = df.Filter('tau_decayMode == {}'.format(dm))
        turnOn_data[dm] = {}
        for wp in working_points:
            wp_bit = ParseEnum(DiscriminatorWP, wp)
            df_wp = df_dm.Filter('({} & (1 << {})) != 0'.format(discr_name, wp_bit))
            turnOn_data[dm][wp] = {}
            for channel in channels:
                turnOn_data[dm][wp][channel] = {}
                df_ch = df_wp.Filter('pass_{} > 0.5'.format(channel))
                for model_name, hist_model in hist_models.items():
                    turn_on = TurnOnData()
                    turn_on.hist_total = df_wp.Histo1D(hist_model, var, branchname_weight)
                    ##print("hist_total:")
                    ##dumpHistogram(turn_on.hist_total)
                    turn_on.hist_passed = df_ch.Histo1D(hist_model, var, branchname_weight)
                    ##print("hist_passed:")
                    ##dumpHistogram(turn_on.hist_passed)
                    turnOn_data[dm][wp][channel][model_name] = turn_on

    return turnOn_data

output_file = ROOT.TFile(args.output + '.root', 'RECREATE')
input_files = [ args.input_data, args.input_dy_mc ]
n_inputs = len(input_files)
branchnames_weight = [ args.branchname_weight_data, args.branchname_weight_dy_mc ] 
labels = [ 'data', 'mc' ]
idx_data = 0
idx_mc   = 1
var = 'tau_pt'
title, x_title = '#tau p_{T}', '#tau p_{T} (GeV)'
decay_modes = args.decay_modes.split(',')
channels = args.channels.split(',')
working_points = args.working_points.split(',')
bins, use_logx = CreateBins(200, False)
bins_fit, _ = CreateBins(200, True)
hist_models = {
    'plot': ROOT.RDF.TH1DModel(var, var, len(bins) - 1, array('d', bins)),
    'fit': ROOT.RDF.TH1DModel(var, var, len(bins_fit) - 1, array('d', bins_fit))
}
turnOn = [ None ] * n_inputs
for input_id in range(n_inputs):
    print("Creating {} histograms...".format(labels[input_id]))
    turnOn[input_id] = CreateHistograms(input_files[input_id], branchnames_weight[input_id], 
                                        channels, decay_modes, 'byDeepTau2017v2p1VSjet',
                                        working_points, hist_models, labels[input_id], var)
dm_labels = {}
for dm in decay_modes:
    if dm == 'all':
        dm_labels[dm] = ''
    else:
        dm_labels[dm] = '_dm{}'.format(dm)
for dm in decay_modes:
    for wp in working_points:
        for channel in channels:
            print('Processing {} {} WP DM = {}'.format(channel, wp, dm))
            for model_name in hist_models.keys():
                turnOn_data = turnOn[idx_data][dm][wp][channel][model_name]
                eff_data = None
                turnOn_mc = turnOn[idx_mc][dm][wp][channel][model_name]
                eff_mc = None
                if 'fit' in model_name:
                    passed_data, total_data, eff_data, passed_mc, total_mc, eff_mc = AutoRebinAndEfficiency(turnOn_data.hist_passed.GetPtr(),
                                                                                                            turnOn_data.hist_total.GetPtr(), 
                                                                                                            turnOn_mc.hist_passed.GetPtr(),
                                                                                                            turnOn_mc.hist_total.GetPtr())
                    #FixEfficiencyBins(passed_data, total_data)
                    #FixEfficiencyBins(passed_mc, total_mc)
                else:
                    passed_data, total_data = turnOn_data.hist_passed.GetPtr(), turnOn_data.hist_total.GetPtr()
                    FixEfficiencyBins(passed_data, total_data)
                    eff_data = ROOT.TEfficiency(passed_data, total_data) 
                    passed_mc, total_mc = turnOn_mc.hist_passed.GetPtr(), turnOn_mc.hist_total.GetPtr()
                    FixEfficiencyBins(passed_mc, total_mc)
                    eff_mc = ROOT.TEfficiency(passed_mc, total_mc)
                name_pattern_data = '{}_{}_{}{}_{}_{{}}'.format(labels[idx_data], channel, wp, dm_labels[dm], model_name)
                turnOn_data.name_pattern = name_pattern_data
                output_file.WriteTObject(total_data, name_pattern_data.format('total'), 'Overwrite')
                output_file.WriteTObject(passed_data, name_pattern_data.format('passed'), 'Overwrite')
                output_file.WriteTObject(eff_data, name_pattern_data.format('eff'), 'Overwrite')
                turnOn_data.eff = eff_data
                name_pattern_mc = '{}_{}_{}{}_{}_{{}}'.format(labels[idx_mc], channel, wp, dm_labels[dm], model_name)
                turnOn_mc.name_pattern = name_pattern_mc
                output_file.WriteTObject(total_mc, name_pattern_mc.format('total'), 'Overwrite')
                output_file.WriteTObject(passed_mc, name_pattern_mc.format('passed'), 'Overwrite')
                output_file.WriteTObject(eff_mc, name_pattern_mc.format('eff'), 'Overwrite')
                turnOn_mc.eff = eff_mc

colors = [ ROOT.kRed, ROOT.kBlack ]
canvas = RootPlotting.CreateCanvas()

n_plots = len(decay_modes) * len(channels) * len(working_points)
plot_id = 0
for channel in channels:
    for wp in working_points:
        for dm in decay_modes:
            if dm == 'all':
                dm_label = ''
                dm_plain_label = ''
            else:
                dm_label = ' DM={}'.format(dm)
                dm_plain_label = '_dm{}'.format(dm)
            ratio_graph = None
            ref_hist = hist_models['plot'].GetHistogram()
            ratio_ref_hist = ref_hist.Clone()
            turnOns = [None] * n_inputs
            curves = [None] * n_inputs
            for input_id in range(n_inputs):
                turnOns[input_id] = turnOn[input_id][dm][wp][channel]['plot']
                curves[input_id] = turnOns[input_id].eff
            y_min, y_max = (0, 1)
            y_title = 'Efficiency'
            title = '{} {}{}'.format(channel, wp, dm_label)
            plain_title = '{}_{}{}'.format(channel, wp, dm_plain_label)
            main_pad, ratio_pad, title_controls = RootPlotting.CreateTwoPadLayout(canvas, ref_hist, ratio_ref_hist,
                                                                                  log_x=use_logx, title=title)
            RootPlotting.ApplyAxisSetup(ref_hist, ratio_ref_hist, x_title=x_title, y_title=y_title,
                                        ratio_y_title='Ratio', y_range=(y_min, y_max * 1.1), max_ratio=1.5)
            legend = RootPlotting.CreateLegend(pos=(0.78, 0.28), size=(0.2, 0.15))
            for input_id in range(n_inputs):
                curve = curves[input_id]
                curve.Draw('SAME')
                RootPlotting.ApplyDefaultLineStyle(curve, colors[input_id])
                legend.AddEntry(curve, labels[input_id], 'PLE')

                if input_id < n_inputs - 1:
                    ratio_graph = RootPlotting.CreateEfficiencyRatioGraph(turnOns[input_id].hist_passed,
                                                                          turnOns[input_id].hist_total,
                                                                          turnOns[-1].hist_passed,
                                                                          turnOns[-1].hist_total)
                    if ratio_graph:
                        output_file.WriteTObject(ratio_graph, 'ratio_{}'.format(plain_title), 'Overwrite')
                        ratio_pad.cd()
                        ratio_color = colors[input_id] if n_inputs > 2 else ROOT.kBlack
                        RootPlotting.ApplyDefaultLineStyle(ratio_graph, ratio_color)
                        ratio_graph.Draw("0PE SAME")
                        main_pad.cd()
            legend.Draw()

            canvas.Update()
            output_file.WriteTObject(canvas, 'canvas_{}'.format(plain_title), 'Overwrite')
            RootPlotting.PrintAndClear(canvas, args.output + '.pdf', plain_title, plot_id, n_plots,
                                       [ main_pad, ratio_pad ])
            plot_id += 1
output_file.Close()
