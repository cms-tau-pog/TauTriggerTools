#!/usr/bin/env python

import argparse
from array import array
import math
import numpy as np
import os
import re
import sys
import ROOT

parser = argparse.ArgumentParser(description='Estimate backgrounds.')
parser.add_argument('--input', required=True, type=str, nargs='+', help="input files")
parser.add_argument('--output-data', required=True, type=str, help="output file prefix for data")
parser.add_argument('--output-dy-mc', required=True, type=str, help="output file prefix for ZTT MC")
parser.add_argument('--output-signal', required=False, type=str, default='signal', help="output file prefix for signal region")
parser.add_argument('--output-w-enriched', required=False, type=str, default='w_enriched', help="output file prefix for w-enriched region")
parser.add_argument('--mode', required=True, type=str, help="subtract backgrounds from data or add backgrounds to ZTT MC")
args = parser.parse_args()

if not(args.mode == "subtract-from-data" or args.mode == "add-to-dy-mc"):
    raise ValueError("Invalid configuration parameter mode = '%s' !!" % args.mode)

path_prefix = '' if 'TauTriggerTools' in os.getcwd() else 'TauTriggerTools/'
sys.path.insert(0, path_prefix + 'Common/python')
from AnalysisTypes import *
from AnalysisTools import *
ROOT.ROOT.EnableImplicitMT(4)
ROOT.gROOT.SetBatch(True)
ROOT.gInterpreter.Declare('#include "{}TauTagAndProbe/interface/PyInterface.h"'.format(path_prefix))

input_vec = ListToStdVector(args.input)

df_input = ROOT.RDataFrame('events', input_vec)
print "df_input = ", df_input
print("sum = %1.2f (%i)" % (df_input.Sum("weight").GetValue(), df_input.Count().GetValue()))
print("")

df_output_data  = None
df_output_dy_mc = None

processes = [ "data", "ztt-mc", "zmm-mc", "w-mc", "ttbar-mc" ]

# define SS->OS extrapolation factor for QCD multijet background 
sf_qcd_SS_to_OS = 1.
print("sf_qcd_SS_to_OS = %1.2f" % sf_qcd_SS_to_OS)
print("")

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

def makeControlPlot(histograms, var, useLogScale, outputFileName):
    
    print("<makeControlPlot>")

    canvasSizeX = 800
    canvasSizeY = 900

    canvas = ROOT.TCanvas("canvas", "", canvasSizeX, canvasSizeY)
    canvas.SetFillColor(10)
    canvas.SetFillStyle(4000)
    canvas.SetFillColor(10)
    canvas.SetTicky()
    canvas.SetBorderSize(2)
    canvas.SetLeftMargin(0.12)
    canvas.SetBottomMargin(0.12)

    topPad = ROOT.TPad("topPad", "topPad", 0.00, 0.35, 1.00, 1.00)
    topPad.SetFillColor(10)
    topPad.SetTopMargin(0.055)
    topPad.SetLeftMargin(0.155)
    topPad.SetBottomMargin(0.030)
    topPad.SetRightMargin(0.050)
    topPad.SetLogy(useLogScale)
  
    bottomPad = ROOT.TPad("bottomPad", "bottomPad", 0.00, 0.00, 1.00, 0.35)
    bottomPad.SetFillColor(10)
    bottomPad.SetTopMargin(0.020)
    bottomPad.SetLeftMargin(0.155)
    bottomPad.SetBottomMargin(0.310)
    bottomPad.SetRightMargin(0.050)
    bottomPad.SetLogy(False)
  
    canvas.cd()
    topPad.Draw()
    topPad.cd()

    histogram_data = histograms['data']
    print("integral['data'] = %1.2f" % histogram_data.Integral())

    xAxis_top = histogram_data.GetXaxis()
    xAxis_top.SetTitle(var);
    xAxis_top.SetTitleOffset(1.2);
    xAxis_top.SetLabelColor(10);
    xAxis_top.SetTitleColor(10);

    yAxis_top = histogram_data.GetYaxis()
    yAxis_top.SetTitle("Events")
    yAxis_top.SetTitleOffset(1.2)
    yAxis_top.SetTitleSize(0.065)
    yAxis_top.SetLabelSize(0.05)
    yAxis_top.SetTickLength(0.04)

    legendTextSize = 0.040
    legendPosX     = 0.740
    legendPosY     = 0.510
    legendSizeX    = 0.190
    legendSizeY    = 0.420

    legend = ROOT.TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC")
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetFillColor(10)
    legend.SetTextSize(legendTextSize)

    histogram_sum = None
    for histogramName, histogram in histograms.items():
        if histogramName != "data":
            if not histogram_sum:
                histogram_sum = histogram.Clone("histogram_sum")
            else:
                histogram_sum.Add(histogram.GetPtr())
    print("integral['sum'] = %1.2f" % histogram_sum.Integral())

    yMin = None
    yMax = None
    if useLogScale:
        yMin = 5.e-1
        yMax = 1.e+1*ROOT.TMath.Max(histogram_data.GetMaximum(), histogram_sum.GetMaximum())
    else:
        yMin = 0.
        yMax = 1.3*ROOT.TMath.Max(histogram_data.GetMaximum(), histogram_sum.GetMaximum())
  
    histogram_data.SetTitle("")
    histogram_data.SetStats(False)
    histogram_data.SetMaximum(yMax)
    histogram_data.SetMinimum(yMin)
    histogram_data.SetMarkerStyle(20)
    markerSize = None
    if histogram_data.GetNbinsX() < 40:
        markerSize = 2
    else:
        markerSize = 1
    histogram_data.SetMarkerSize(markerSize)
    histogram_data.SetMarkerColor(1)
    histogram_data.SetLineColor(1)
    legend.AddEntry(histogram_data.GetPtr(), "observed", "p");
    histogram_data.Draw("ep");

    colors = {}
    colors['ztt-mc']   = 796
    colors['zmm-mc']   = 842
    colors['w-mc']     = 634
    colors['ttbar-mc'] = 592
    colors['qcd']      = 606

    legendEntries = {}
    legendEntries['ztt-mc']   = "Z#rightarrow#tau#tau"
    legendEntries['zmm-mc']   = "Z#rightarrow#mu#mu"
    legendEntries['w-mc']     = "W+jets"
    legendEntries['ttbar-mc'] = "t#bar{t}+jets"
    legendEntries['qcd']      = "Multijet"

    histograms_stack = ROOT.THStack("stack", "");
    for histogramName in [ "qcd", "w-mc", "ttbar-mc", "zmm-mc", "ztt-mc" ]:
        histogram = histograms[histogramName]
        makeBinContentsPositive(histogramName, histogram)
        print("integral['%s'] = %1.2f" % (histogramName, histogram.Integral()))
        histogram.SetFillColor(colors[histogramName])
        histogram.SetLineColor(1)
        histograms_stack.Add(histogram.GetPtr())
    for histogramName in reversed([ "qcd", "w-mc", "ttbar-mc", "zmm-mc", "ztt-mc" ]):
        histogram = histograms[histogramName]
        legend.AddEntry(histogram.GetPtr(), legendEntries[histogramName], "f")
    histograms_stack.Draw("histsame")

    histogram_data.Draw("epsame")
    histogram_data.Draw("axissame")
   
    legend.Draw()
    
    canvas.cd()
    bottomPad.Draw()
    bottomPad.cd()
  
    histogram_ratio = histogram_data.Clone("histogram_ratio")
    histogram_ratio.Reset()
    if not histogram_ratio.GetSumw2N():
        histogram_ratio.Sumw2()
    histogram_ratio.Divide(histogram_data.GetPtr(), histogram_sum);
    for i in range(histogram_ratio.GetNbinsX()):
        binContent = histogram_ratio.GetBinContent(i + 1)
        histogram_ratio.SetBinContent(i + 1, binContent - 1.0)
    histogram_ratio.SetTitle("")
    histogram_ratio.SetStats(False)
    histogram_ratio.SetMinimum(-0.50)
    histogram_ratio.SetMaximum(+0.50)
    histogram_ratio.SetMarkerStyle(histogram_data.GetMarkerStyle())
    histogram_ratio.SetMarkerSize(histogram_data.GetMarkerSize())
    histogram_ratio.SetMarkerColor(histogram_data.GetMarkerColor())
    histogram_ratio.SetLineColor(histogram_data.GetLineColor())
    
    xAxis_bottom = histogram_ratio.GetXaxis()
    xAxis_bottom.SetTitle(xAxis_top.GetTitle())
    xAxis_bottom.SetLabelColor(1)
    xAxis_bottom.SetTitleColor(1)
    xAxis_bottom.SetTitleOffset(1.20)
    xAxis_bottom.SetTitleSize(0.12)
    xAxis_bottom.SetLabelOffset(0.02)
    xAxis_bottom.SetLabelSize(0.10)
    xAxis_bottom.SetTickLength(0.055)
    
    yAxis_bottom = histogram_ratio.GetYaxis()
    yAxis_bottom.SetTitle("#frac{Data - Expectation}{Expectation}")
    yAxis_bottom.SetTitleOffset(0.80)
    yAxis_bottom.SetNdivisions(505)
    yAxis_bottom.CenterTitle()
    yAxis_bottom.SetTitleSize(0.09)
    yAxis_bottom.SetLabelSize(0.10)
    yAxis_bottom.SetTickLength(0.04)
    
    histogram_ratio.Draw("ep")
    
    line = ROOT.TF1("line","0", xAxis_bottom.GetXmin(), xAxis_bottom.GetXmax())
    line.SetLineStyle(3)
    line.SetLineWidth(1)
    line.SetLineColor(1)
    line.Draw("same")
        
    histogram_ratio.Draw("epsame")
  
    canvas.Update()
    idx = outputFileName.rfind('.')
    outputFileName_plot = outputFileName[0:idx]
    if useLogScale: 
        outputFileName_plot += "_log"
    else:
        outputFileName_plot += "_linear"
    canvas.Print(outputFileName_plot + ".png")
    canvas.Print(outputFileName_plot + ".pdf")
#----------------------------------------------------------------------------------------------------

#-------------- MAKE ROOT FILES FOR SIGNAL AND W-ENRICHED SIDEBANDS (only in "subtract-from-data" mode)
if args.mode == "subtract-from-data":
    df_signal_data = df_input.Filter("selection == {} && type == {}".format(selection_signal, type_data))
    df_signal_ztt_mc = df_input.Filter("selection == {} && type == {}".format(selection_signal, type_ztt_mc))
    df_signal_qcd = df_input.Filter("selection == {}".format(selection_SS_low_mT))
    df_w_enriched_data = df_input.Filter("selection == {} && type == {}".format(selection_w_enriched, type_data))
    print("Total weight: df_w_enriched_data = ", df_w_enriched_data.Sum("weight").GetValue())
    df_w_enriched_mc = df_input.Filter('''
                                       (selection == {} && type == {}) 
                                       || (selection == {} && type == {}) 
                                       || (selection == {} && type == {})
                                       '''.format(selection_w_enriched, type_ttbar_mc, 
                                                  selection_w_enriched, type_w_mc, 
                                                  selection_w_enriched, type_zmm_mc))
    print("Total weight: df_w_enriched_mc = ", df_w_enriched_mc.Sum("weight").GetValue())
    df_signal_data.Snapshot('events', args.output_signal + '_data.root')
    df_signal_ztt_mc.Snapshot('events', args.output_signal + '_ztt_mc.root')
    df_signal_qcd.Snapshot('events', args.output_signal + '_signal_qcd_inputs.root')
    df_w_enriched_data.Snapshot('events', args.output_w_enriched + '_data.root')
    df_w_enriched_mc.Snapshot('events', args.output_w_enriched + '_mc.root')
#----------------------------------------------------------------------------------------------------

# step 1: determine scale-factor for W+jets background in SS region
df_SS_high_mT = df_input.Filter("selection == %i" % selection_SS_high_mT)
sum_SS_high_mT = {}
for process in processes:
    df_SS_high_mT_process = df_SS_high_mT.Filter("type == %i" % get_type(process))
    sum_SS_high_mT[process] = df_SS_high_mT_process.Sum("weight").GetValue()
    print("sum_SS_high_mT['%s'] = %1.2f (%i)" % (process, sum_SS_high_mT[process], df_SS_high_mT_process.Count().GetValue()))
sf_w_mc_SS = (sum_SS_high_mT['data'] - (sum_SS_high_mT["ztt-mc"] + sum_SS_high_mT['zmm-mc'] + sum_SS_high_mT['ttbar-mc']))/sum_SS_high_mT['w-mc']
print("sf_w_mc_SS = %1.2f" % sf_w_mc_SS)
print("")

# step 2: determine QCD multijet background in SS region 
#        (Note: QCD multijet background in SS high mT sideband assumed to be negligible)
df_SS_low_mT = df_input.Filter("selection == %i" % selection_SS_low_mT)
sum_SS_low_mT = {}
sum_data_SS_low_mT = 0.
sum_mc_SS_low_mT = 0.
for process in processes:
    df_SS_low_mT_process = df_SS_low_mT.Filter("type == %i" % get_type(process))
    sum_SS_low_mT[process] = df_SS_low_mT_process.Sum("weight").GetValue()
    if process == "w-mc":
        sum_SS_low_mT[process] *= sf_w_mc_SS
    print("sum_SS_low_mT['%s'] = %1.2f (%i)" % (process, sum_SS_low_mT[process], df_SS_low_mT_process.Count().GetValue()))
    if process == "data":
        sum_data_SS_low_mT += sum_SS_low_mT[process]
    else:
        sum_mc_SS_low_mT += sum_SS_low_mT[process]
sum_SS_low_mT['qcd'] = sum_data_SS_low_mT - sum_mc_SS_low_mT
print("sum_SS_low_mT['qcd'] = %1.2f" % sum_SS_low_mT['qcd'])
print("")

# step 3: determine scale-factor for W+jets background in OS region 
#        (Note: QCD multijet background in OS high mT sideband assumed to be negligible)
df_OS_high_mT = df_input.Filter("selection == %i" % selection_OS_high_mT)
sum_OS_high_mT = {}
for process in processes:
    df_OS_high_mT_process = df_OS_high_mT.Filter("type == %i" % get_type(process))
    sum_OS_high_mT[process] = df_OS_high_mT_process.Sum("weight").GetValue()
    print("sum_OS_high_mT['%s'] = %1.2f (%i)" % (process, sum_OS_high_mT[process], df_OS_high_mT_process.Count().GetValue()))
sf_w_mc_OS = (sum_OS_high_mT['data'] - (sum_OS_high_mT["ztt-mc"] + sum_OS_high_mT['zmm-mc'] + sum_OS_high_mT['ttbar-mc']))/sum_OS_high_mT['w-mc']
print("sf_w_mc_OS = %1.2f" % sf_w_mc_OS)
print("")

# step 4: print event yields in "signal" region for input RDataFrame objects
df_OS_low_mT = df_input.Filter("selection == %i" % selection_OS_low_mT)
sum_OS_low_mT = {}
sum_data_OS_low_mT = 0.
sum_mc_OS_low_mT = 0.
for process in processes:
    df_OS_low_mT_process = df_OS_low_mT.Filter("type == %i" % get_type(process))
    sum_OS_low_mT[process] = df_OS_low_mT_process.Sum("weight").GetValue()
    if process == "w-mc":
        sum_OS_low_mT[process] *= sf_w_mc_OS
    print("sum_OS_low_mT['%s'] = %1.2f (%i)" % (process, sum_OS_low_mT[process], df_OS_low_mT_process.Count().GetValue()))
    if process == "data":
        sum_data_OS_low_mT += sum_OS_low_mT[process]
    else:
        sum_mc_OS_low_mT += sum_OS_low_mT[process]
print("sum_OS_low_mT['qcd'] = %1.2f" % (sum_SS_low_mT['qcd']*sf_qcd_SS_to_OS))
print("")

# step 5: build RDataFrame object for 'data'
final_weight_data = ROOT.final_weight_data.Initialize(sf_qcd_SS_to_OS, sf_w_mc_OS, sf_w_mc_SS)
if args.mode == "subtract-from-data":
    df_output_data = df_input.Filter("(selection == %i && type == %i) || (selection == %i) || (selection == %i && (type == %i || type == %i || type == %i))" % (selection_OS_low_mT, type_data, selection_SS_low_mT, selection_OS_low_mT, type_zmm_mc, type_w_mc, type_ttbar_mc))
    df_output_data = df_output_data.Define("final_weight", "final_weight_data::GetDefault().operator()(selection, type, weight)")
elif args.mode == "add-to-dy-mc":
    df_output_data = df_input.Filter("selection == %i && type == %i" % (selection_OS_low_mT, type_data))

# step 6: build RDataFrame object for 'dy-mc'
final_weight_dy_mc = ROOT.final_weight_dy_mc.Initialize(sf_qcd_SS_to_OS, sf_w_mc_OS, sf_w_mc_SS)
if args.mode == "subtract-from-data":
    df_output_dy_mc = df_input.Filter("selection == %i && type == %i" % (selection_OS_low_mT, type_ztt_mc))
elif args.mode == "add-to-dy-mc":
    df_output_dy_mc = df_input.Filter("(selection == %i && type == %i) || (selection == %i) || (selection == %i && (type == %i || type == %i || type == %i))" % (selection_OS_low_mT, type_ztt_mc, selection_SS_low_mT, selection_OS_low_mT, type_zmm_mc, type_w_mc, type_ttbar_mc))
    df_output_dy_mc = df_output_dy_mc.Define("final_weight", "final_weight_dy_mc::GetDefault().operator()(selection, type, weight)")

# step 7: print data and dy-mc event yields in "signal" region for output RDataFrame objects
print "df_output_data = ", df_output_data
weight_data = None
if args.mode == "subtract-from-data":
    weight_data = "final_weight"
else:
    weight_data = "weight"
sum_OS_low_mT_data = df_output_data.Sum(weight_data).GetValue()
print("sum_OS_low_mT_data = %1.2f (%i)" % (sum_OS_low_mT_data, df_output_data.Count().GetValue()))
print "df_output_dy_mc = ", df_output_dy_mc
weight_dy_mc = None
if args.mode == "add-to-dy-mc":
    weight_dy_mc = "final_weight"
else:
    weight_dy_mc = "weight"
sum_OS_low_mT_dy_mc = df_output_dy_mc.Sum(weight_dy_mc).GetValue()
print("sum_OS_low_mT_dy_mc = %1.2f (%i)" % (sum_OS_low_mT_dy_mc, df_output_dy_mc.Count().GetValue()))
print("")

# step 8: make control plots (only implemented for "add-to-dy-mc" mode so far)
var = "tau_pt"
branchname_weight_data  = "weight"
branchname_weight_dy_mc = "final_weight"
hist_model = ROOT.RDF.TH1DModel(var, var, 18, 20., 200.)
histograms = {}
if args.mode == "add-to-dy-mc":
    discr_name = "byDeepTau2017v2p1VSjet"
    for wp in [ "VVVLoose", "VVLoose", "VLoose", "Loose", "Medium", "Tight", "VTight", "VVTight" ]:
        wp_bit = ParseEnum(DiscriminatorWP, wp)
        offlineTauSel  = "tau_pt > 20 && abs(tau_eta) < 2.3"
        offlineTauSel += " && (tau_decayMode == 0 || tau_decayMode == 1 || tau_decayMode == 2 || tau_decayMode == 10 || tau_decayMode == 11)"
        offlineTauSel += " && (%s & (1 << %i))" % (discr_name, wp_bit)
        df_data_passing_offlineTauSel = df_output_data.Filter(offlineTauSel)
        histograms['data']     = df_data_passing_offlineTauSel.Histo1D(hist_model, var, branchname_weight_data)
        df_dy_mc_passing_offlineTauSel = df_output_dy_mc.Filter(offlineTauSel)
        histograms['ztt-mc']   = df_dy_mc_passing_offlineTauSel.Filter("selection == %i && type == %i" % (selection_OS_low_mT, type_ztt_mc)).Histo1D(hist_model, var, branchname_weight_dy_mc)
        histograms['zmm-mc']   = df_dy_mc_passing_offlineTauSel.Filter("selection == %i && type == %i" % (selection_OS_low_mT, type_zmm_mc)).Histo1D(hist_model, var, branchname_weight_dy_mc)
        histograms['w-mc']     = df_dy_mc_passing_offlineTauSel.Filter("selection == %i && type == %i" % (selection_OS_low_mT, type_w_mc)).Histo1D(hist_model, var, branchname_weight_dy_mc)
        histograms['ttbar-mc'] = df_dy_mc_passing_offlineTauSel.Filter("selection == %i && type == %i" % (selection_OS_low_mT, type_ttbar_mc)).Histo1D(hist_model, var, branchname_weight_dy_mc)
        histograms['qcd']      = df_dy_mc_passing_offlineTauSel.Filter("selection == %i" % selection_SS_low_mT).Histo1D(hist_model, var, branchname_weight_dy_mc)
        #outputFileName = "estimateBackgrounds_%s%s_%s.pdf" % (discr_name, wp, var)
        outputFileName = "%s_%s%s_%s.pdf" % (args.output_dy_mc, discr_name, wp, var)
        makeControlPlot(histograms, var, True, outputFileName)
        print("")

# step 9: write RDataFrame objects to output files
df_output_data.Snapshot('events', args.output_data + '.root')
df_output_dy_mc.Snapshot('events', args.output_dy_mc + '.root')
