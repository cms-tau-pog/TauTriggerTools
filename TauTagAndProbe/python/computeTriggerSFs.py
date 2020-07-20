#!/usr/bin/env python                                                                                                                                                                                   

import argparse
from array import array
import math
import numpy as np
import os
import re
import sys
import ROOT

parser = argparse.ArgumentParser(description='Estimate QCD backgrounds.')
parser.add_argument('--input_stage2p5', required=True, type=str, default='', help="input Stage 2.5 file")
parser.add_argument('--input_signal', required=True, type=str, default='', help="input Signal region root file with fitted turn ons")
parser.add_argument('--input_w_enriched', required=True, type=str, default='', help="input W-enriched region root file w fitted turn ons")
parser.add_argument('--channels', required=False, type=str, default='etau,mutau,ditau', help="channels to process")
parser.add_argument('--decay-modes', required=False, type=str, default='all,0,1,10,11', help="decay modes to process")
parser.add_argument('--working-points', required=False, type=str,
                    default='VVVLoose,VVLoose,VLoose,Loose,Medium,Tight,VTight,VVTight',
                    help="working points to process")
parser.add_argument('--output', required=True, type=str, help="output file name")
args = parser.parse_args()

if ((args.input_stage2p5 == '') or (args.input_signal == '') or (args.input_w_enriched == '')):
    raise ValueError("One or more of these configuration parameters are Invalid = '%s', '%s', '%s' !!" % (args.input_stage2p5, args.input_signal, args.input_w_enriched))

path_prefix = '' if 'TauTriggerTools' in os.getcwd() else 'TauTriggerTools/'
sys.path.insert(0, path_prefix + 'Common/python')
from AnalysisTypes import *
from AnalysisTools import *
import RootPlotting
ROOT.ROOT.EnableImplicitMT(4)
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2() # CV: This does not seem to work... all bin-error are the exact square-root of the bin-contents !!   
ROOT.gInterpreter.Declare('#include "{}TauTagAndProbe/interface/PyInterface.h"'.format(path_prefix))
RootPlotting.ApplyDefaultGlobalStyle()

#input_vec = ListToStdVector(args.input_stage2p5)
#processes = [ "data", "ztt-mc", "zmm-mc", "w-mc", "ttbar-mc" ]
#qcd_input_processes = [ "data", "zmm-mc", "w-mc", "ttbar-mc" ]

#----------------------------------------------------------------------------------------------------                                                                           
# define integer constants                                                                                                                                                      
type_data            = 0
type_ztt_mc          = 1
type_zmm_mc          = 2
type_w_mc            = 3
type_ttbar_mc        = 4

selection_OS_low_mT  = 0
selection_OS_high_mT = 1
selection_SS_low_mT  = 2
selection_SS_high_mT = 3
selection_signal = 4
selection_w_enriched = 5
#----------------------------------------------------------------------------------------------------

def get_type(process):
    if process == "data":
        return type_data
    elif process == "ztt-mc":
        return type_ztt_mc
    elif process == "zmm-mc":
        return type_zmm_mc
    elif process == "w-mc":
        return type_w_mc
    elif process == "ttbar-mc":
        return type_ttbar_mc
    else:
        raise ValueError("Invalid function argument: process = '%s' !!" % process)
#----------------------------------------------------------------------------------------------------                                                                           

def makeBinContentsPositive(histogramName, histogram):
    print("<makeBinContentsPositive>: histogram = %s" % histogramName)
    integral_original = histogram.Integral()
    for i in range(histogram.GetNbinsX()):
        binContent = histogram.GetBinContent(i + 1)
        binError = histogram.GetBinError(i + 1)
        if binContent < 0.:
            histogram.SetBinContent(i + 1, 0.)
            histogram.SetBinError(i + 1, math.sqrt(binContent**2 + binError**2))
    integral_modified = histogram.Integral()
    print("integral: original = %1.2f, modified = %1.2f" % (integral_original, integral_modified))
    if integral_modified > 0.:
        histogram.Scale(integral_original/integral_modified)

def extractBinEdges(histogramName, histogram):
    print("<extractBinEdges>: histogram = %s" % histogramName)
    print(" histogram.GetNbinsX() ", histogram.GetNbinsX())
    #bin_array = np.empty(histogram.GetNbinsX()) ## Creating an empty numpy array of length "histogram.GetNbinsX()"
    #bin_array = np.zeros((histogram.GetNbinsX(),)) ## Creating numpy array full of zeros of length "histogram.GetNbinsX()"
    bin_list = [0.]
    for i in range(histogram.GetNbinsX()):
        hist_BinLowEdge = histogram.GetBinLowEdge(i + 1)
        #bin_array.append(hist_BinLowEdge)
        #bin_array[i+1] = hist_BinLowEdge
        bin_list.append(hist_BinLowEdge)
    bin_array = np.array(bin_list)    
    print("bin_array ", bin_array)    
    return bin_array    


class AlphaInfo:
    def __init__(self):
        self.hist_alpha = None
        self.hist_one_minus_alpha = None
        self.eff = None

def IncreaseBinning(alpha_orig, AlphaInfo_Obj, var, bins):
    alpha = ROOT.TH1F("alpha", var, len(bins) - 1, array('d', bins)) 
    one_minus_alpha = alpha.Clone()
    for n in range(1, alpha.GetNbinsX() + 1):
        x = alpha.GetBinCenter(n)
        bin_orig = alpha_orig.FindFixBin(x)
        alpha.SetBinContent(n, alpha_orig.GetBinContent(bin_orig))
        alpha.SetBinError(n, alpha_orig.GetBinError(bin_orig))
        one_minus_alpha.SetBinContent(n, 1 - alpha_orig.GetBinContent(bin_orig))
        one_minus_alpha.SetBinError(n, alpha_orig.GetBinError(bin_orig))
    AlphaInfo_Obj.hist_alpha = alpha
    AlphaInfo_Obj.hist_one_minus_alpha = one_minus_alpha

def ComputeQCD(df, hist_model, var, branchname_weight, sf_qcd_SS_to_OS = 1.0):
    df_SS = df.Filter("selection == {}".format(selection_SS_low_mT))
    #hist_data = (df_SS.Filter("type == {}".format(type_data))).Histo1D(hist_model, var, branchname_weight)
    hist_qcd = (df_SS.Filter("type == {}".format(type_data))).Histo1D(hist_model, var, branchname_weight)
    hist_mc = (df_SS.Filter('''
                            type == {} or type == {} or type == {} or type == {}
                            '''.format(type_ztt_mc, type_zmm_mc, type_w_mc, type_ttbar_mc))).Histo1D(hist_model, var, branchname_weight)
    #hist_qcd = hist_data.Clone()
    print("SS Data :", hist_qcd.Integral())
    print("SS MC :", hist_mc.Integral())
    #print("type(hist_qcd) :", type(hist_qcd))
    #print("type(hist_mc) :", type(hist_mc))
    hist_qcd.Add(hist_mc.GetPtr(), -1)
    hist_qcd.Scale(sf_qcd_SS_to_OS)
    print("QCD yield computed :", hist_qcd.Integral())
    if(hist_qcd.Integral() < 0.):
        print("Data driven QCD yield is -ve, setting it to zero by hand")
        hist_qcd.Reset()
    #makeBinContentsPositive("hist_qcd", hist_qcd)
    print("QCD final yield :", hist_qcd.Integral())
    return hist_qcd

def Compute_eff_data_true(AlphaInfo_Obj, eff_data_sgn, eff_data_fake):
    eff_data_true_tmp = eff_data_fake.Clone()
    eff_data_true = eff_data_fake.Clone()
    eff_data_true.Reset()
    one_minus_alpha = AlphaInfo_Obj.hist_one_minus_alpha.Clone()
    alpha = AlphaInfo_Obj.hist_alpha.Clone()
    print("type(eff_data_true_tmp) ", type(eff_data_true_tmp))
    print("type(one_minus_alpha) ", type(one_minus_alpha))
    print("type(alpha) ", type(alpha))
    eff_data_true.Multiply(eff_data_true_tmp, one_minus_alpha, 1, -1) 
    #eff_data_true.Multiply(one_minus_alpha, -1)
    eff_data_true.Add(eff_data_sgn)
    eff_data_true.Divide(alpha)
    eff_data_true.SetName("eff_data_true")
    return eff_data_true

def ComputeSF(input_files, output_file, branchname_weight,
                 channels, decay_modes, discr_name,
                 working_points, hist_models, label, 
                 var, sf_qcd_SS_to_OS):
    print("<ComputeAlpha>:")
    print("Stage 2.5 input_file = '%s'" % input_files[0]) 
    df = ROOT.RDataFrame('events', input_files[0])
    print("Fitted TurnOn signal region input file = '%s'" % input_files[1]) 
    file_signal = ROOT.TFile(input_files[1], "READ")
    print("Fitted TurnOn w-enriched region input file = '%s'" % input_files[2]) 
    file_w_enriched = ROOT.TFile(input_files[2], "READ")
    #alpha = {}
    dm_labels = {}
    turnon_label = ""
    for dm in decay_modes:
        print("decay mode :", dm)
        if dm == 'all':
            dm_labels[dm] = ''
            df_dm = df
        else:
            dm_labels[dm] = '_dm{}'.format(dm)
            df_dm = df.Filter('tau_decayMode == {}'.format(dm))
        #alpha[dm] = {}
        for wp in working_points:
            print("Working point :", wp)
            wp_bit = ParseEnum(DiscriminatorWP, wp)
            df_wp = df_dm.Filter('({} & (1 << {})) != 0'.format(discr_name, wp_bit))
            #alpha[dm][wp] = {}
            for channel in channels:
                print("channel :", channel)
                turnon_label = "{}_{}_dm{}_fitted".format(channel, wp, dm)
                print("turnon_label: ", turnon_label)
                #alpha[dm][wp][channel] = {}
                df_ch = df_wp.Filter('pass_{} > 0.5'.format(channel))
                for model_name, hist_model in hist_models.items():
                    A = AlphaInfo()
                    hist_qcd = ComputeQCD(df_ch, hist_model, var, branchname_weight, sf_qcd_SS_to_OS) ## TO BE IMPLEMENTED
                    hist_signal = (df_ch.Filter('''
                                                selection == {} && type == {}
                                                '''.format(selection_signal, type_ztt_mc))).Histo1D(hist_model, var, branchname_weight)
                    hist_mc_bkg = (df_ch.Filter('''
                                                selection == {} && ((type == {}) or (type == {}) or (type == {}))
                                                '''.format(selection_signal, type_w_mc, type_ttbar_mc, type_zmm_mc))).Histo1D(hist_model, var, branchname_weight)
                    hist_total = hist_mc_bkg.Clone() ## W_mc + TT_mc + Zmm_mc
                    print("W_mc + TT_mc + Zmm_mc : ", hist_total.Integral())
                    #print("type(hist_qcd) :", type(hist_qcd))
                    #print("type(hist_total) :", type(hist_total))
                    hist_total.Add(hist_qcd.GetPtr()) ## Total Bg = W_mc + TT_mc + Zmm_mc + QCD
                    print("W_mc + TT_mc + Zmm_mc + QCD : ", hist_total.Integral())
                    alpha_orig = hist_signal.Clone() ## Numerator for alpha_orig = Signal = Ztt_mc
                    print("Ztt_mc : ", hist_signal.Integral())
                    hist_total.Add(hist_signal.GetPtr()) ## Denominator for alpha_orig = Signal + Total Bg
                    print("(W_mc + TT_mc + Zmm_mc + QCD) + Ztt_mc : ", hist_total.Integral())
                    #print("type(alpha_orig) :", type(alpha_orig))
                    #print("type(hist_total) :", type(hist_total))
                    alpha_orig.Divide(hist_total) ## alpha_orig = Signal/(Signal + Total Bg)
                    #hist_signal.Divide(hist_total.GetPtr()) ## Signal/(Signal + Total Bg)

                    mc_turnon_label = "mc_" + turnon_label
                    data_turnon_label = "data_" + turnon_label

                    eff_mc_true = file_signal.Get(mc_turnon_label)
                    eff_mc_true.SetName("eff_mc_true")
                    eff_mc_true_orig = eff_mc_true.Clone()

                    eff_data_sgn = file_signal.Get(data_turnon_label)
                    eff_data_sgn.SetName("eff_data_sgn")
                    eff_data_sgn_orig = eff_data_sgn.Clone() 

                    eff_mc_fake = file_w_enriched.Get(mc_turnon_label)
                    print("type(eff_mc_fake) ", type(eff_mc_fake))
                    print("eff_mc_fake.Integral() ", eff_mc_fake.Integral())
                    eff_mc_fake.SetName("eff_mc_fake")
                    eff_mc_fake_orig = eff_mc_fake.Clone()  

                    eff_data_fake = file_w_enriched.Get(data_turnon_label)
                    eff_data_fake.SetName("eff_data_fake")
                    eff_data_fake_orig = eff_data_fake.Clone()
                    #sf_fake = eff_data_fake.Clone()
                    #sf_fake.Reset()

                    print("type(eff_mc_true) ", type(eff_mc_true))
                    print("eff_mc_true.Integral() ", eff_mc_true.Integral())
                    print("type(eff_mc_fake) ", type(eff_mc_fake))
                    print("eff_mc_fake.Integral() ", eff_mc_fake.Integral())

                    bins_w_enriched = extractBinEdges(eff_data_fake.GetName(), eff_data_fake)
                    bins_signal = extractBinEdges(eff_data_sgn.GetName(), eff_data_sgn)
                    if(np.array_equal(bins_signal,bins_w_enriched)):
                        print("Signal and W-enriched arrays are identical")
                        bins = bins_w_enriched
                        print("bins_signal ", bins_signal)
                    else:    
                        print("Signal and W-enriched arrays are not equal")
                        bins = np.arange(20., 1000., 0.1) ## 19800 Bins of 0.1 GeV spacing from 20-1000 GeV

                    IncreaseBinning(alpha_orig, A, var, bins)
                    eff_data_true = Compute_eff_data_true(A, eff_data_sgn, eff_data_fake)
                    eff_data_true_orig = eff_data_true.Clone()
                    #sf_true = eff_data_true.Clone()
                    #sf_true.Reset()
                    eff_data_true.Divide(eff_mc_true)
                    eff_data_true.SetName("sf_true_taus")
                    eff_data_fake.Divide(eff_mc_fake)
                    eff_data_fake.SetName("sf_fake_taus")
                    out_name_pattern = '{}_{{}}'.format(turnon_label)
                    output_file.WriteTObject(eff_data_true, out_name_pattern.format('sf_true'), 'Overwrite')
                    output_file.WriteTObject(eff_data_fake, out_name_pattern.format('sf_fake'), 'Overwrite')
                    output_file.WriteTObject(A.hist_alpha, out_name_pattern.format('alpha'), 'Overwrite')
                    output_file.WriteTObject(A.hist_one_minus_alpha, out_name_pattern.format('one_minus_alpha'), 'Overwrite')
                    output_file.WriteTObject(eff_data_sgn_orig, out_name_pattern.format('eff_data_sgn'), 'Overwrite')
                    output_file.WriteTObject(eff_data_true_orig, out_name_pattern.format('eff_data_true'), 'Overwrite')
                    output_file.WriteTObject(eff_data_fake_orig, out_name_pattern.format('eff_data_fake'), 'Overwrite')
                    output_file.WriteTObject(eff_mc_true_orig, out_name_pattern.format('eff_mc_true'), 'Overwrite')
                    output_file.WriteTObject(eff_mc_fake_orig, out_name_pattern.format('eff_mc_fake'), 'Overwrite')
                    #alpha[dm][wp][channel][model_name] = A
    file_signal.Close() 
    file_w_enriched.Close()
    output_file.Close()
    print('All done.')
    #return alpha
                

#input_file = ROOT.TFile(args.input_stage2p5, 'READ')
#output_file = ROOT.TFile(args.output + '.root', 'RECREATE')
output_file = ROOT.TFile('{}.root'.format(args.output), 'RECREATE', '', ROOT.RCompressionSetting.EDefaults.kUseSmallest)
#df_input = ROOT.RDataFrame('events', input_vec)
input_files = [args.input_stage2p5, args.input_signal, args.input_w_enriched]
n_inputs = len(input_files) - 2
branchnames_weight = ['weight']
#labels = [ 'data', 'mc', 'qcd' ]
#idx_data = 0
#idx_mc   = 1
#idx_qcd = 2
labels = ['alpha' ]
idx_alpha = 0
sf_qcd_SS_to_OS = 1.0
var = 'tau_pt'
title, x_title = '#tau p_{T}', '#tau p_{T} (GeV)'
channels = args.channels.split(',')
decay_modes = args.decay_modes.split(',')
working_points = args.working_points.split(',')
bins = np.arange(20., 1000., 0.1) ## 1800 Bins of 0.1 GeV spacing from 20-1000 GeV
print("type(bins): ", type(bins))
hist_models = {
    'plot': ROOT.RDF.TH1DModel(var, var, len(bins) - 1, array('d', bins))
}


alpha = [ None ] * n_inputs  ## List of dictionaries
for input_id in range(n_inputs):
    print("Creating {} histograms...".format(labels[input_id]))
    alpha[input_id] = ComputeSF(input_files, output_file, branchnames_weight[input_id],
                                   channels, decay_modes, 'byDeepTau2017v2p1VSjet',
                                   working_points, hist_models, labels[input_id], 
                                   var, sf_qcd_SS_to_OS)
