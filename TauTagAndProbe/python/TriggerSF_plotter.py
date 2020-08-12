import argparse
from array import array
import math
import numpy as np
import os
import re
import sys
import ROOT
from ROOT import TAttFill

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)

# python TriggerSF_plotter.py --era 2016 --decay-mode 0 --channels ditau --working-points Medium  --inputFilePath-new $PWD --outputFilePath $PWD/Tau_Trigger_sf_plots

parser = argparse.ArgumentParser(description='Plotter for Trigger SFs')
#parser.add_argument('--inputFilePath-old', required=True, type=str, help="input file Path for the Konstantin's old files")
parser.add_argument('--inputFilePath-new', required=True, type=str, help="input file Path for the files")
parser.add_argument('--outputFilePath', required=True, type=str, help="Name of the output file path")
parser.add_argument('--era', required=False, type=str, default='2016,2017,2018', help="Era")
parser.add_argument('--decay-mode', required=False, type=str, default='all,0,1,10,11', help="decay mode indices")
parser.add_argument('--channels', required=False, type=str, default='etau,mutau,ditau', help="channels to process")
parser.add_argument('--working-points', required=False, type=str,
                    default='VVVLoose,VVLoose,VLoose,Loose,Medium,Tight,VTight,VVTight',
                    help="working points to process")
args = parser.parse_args()

#InputFilePath_old = args.inputFilePath_old
InputFilePath_new = args.inputFilePath_new
OutputFilePath = args.outputFilePath
Eras = args.era.split(',')
Decay_modes = args.decay_mode.split(',')
Channels = args.channels.split(',')
Working_points = args.working_points.split(',')


def makePlot_NewSFs(histo_dict, outputFileName):
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
    canvas.cd() 

    histogram_true = histo_dict['sf_true']
    histogram_fake = histo_dict['sf_fake']
    
    histogram_true.SetMaximum(2.0)
    histogram_fake.SetMaximum(2.0)
    
    histogram_true.SetLineColor(2)
    histogram_fake.SetLineColor(3)

    histogram_true.SetMarkerColor(2)
    histogram_fake.SetMarkerColor(3)
    
    histogram_true.SetMarkerStyle(20)
    histogram_fake.SetMarkerStyle(21)
    
    histogram_true.SetMarkerSize(0.7)
    histogram_fake.SetMarkerSize(0.7)
    
    histogram_true.SetFillColor(2)
    histogram_fake.SetFillColor(3)
    
    histogram_true.SetFillStyle(3004)
    histogram_fake.SetFillStyle(3005)


    xAxis_top = histogram_true.GetXaxis()
    xAxis_top.SetTitle("Tau p_{T} (GeV)");
    xAxis_top.SetTitleOffset(1.2);
    xAxis_top.SetTitleSize(0.03)
    xAxis_top.SetLabelSize(0.03)


    yAxis_top = histogram_true.GetYaxis()
    yAxis_top.SetTitle("Data/MC SF")
    yAxis_top.SetTitleOffset(1.2)
    yAxis_top.SetTitleSize(0.03)
    yAxis_top.SetLabelSize(0.03)
    yAxis_top.SetTickLength(0.04)
    
    legendTextSize = 0.040
    legendPosX     = 0.740
    legendPosY     = 0.510
    legendSizeX    = 0.190
    legendSizeY    = 0.420

    legend = ROOT.TLegend(0.7, 0.7, 0.85, 0.85, "", "brNDC")
    legend.SetFillStyle(0)
    legend.SetFillColor(10)
    legend.SetTextSize(0.018)

    legend.AddEntry(histogram_true, "True taus", "f")
    legend.AddEntry(histogram_fake, "Fake taus", "f")

    histogram_true.GetXaxis().SetRangeUser(20.,200.)
    histogram_fake.GetXaxis().SetRangeUser(20.,200.)

    histogram_true.GetYaxis().SetRangeUser(0., 1.4)
    histogram_fake.GetYaxis().SetRangeUser(0., 1.4)    

    histogram_true.Draw('E2') 
    histogram_fake.Draw('E2 same') 
    
    legend.Draw()

    canvas.Update()
    canvas.Print(outputFileName + ".pdf")
    canvas.Print(outputFileName + ".png")
    canvas.Print(outputFileName + ".root")


def makePlot_OldSFs(histo_dict, outputFileName):
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

    canvas.cd() 

    #histogram_old = histo_dict['old']
    histogram_new1 = histo_dict['sub_from_data']
    histogram_new2 = histo_dict['add_to_dy_mc']

    #histogram_old.SetMaximum(2.0)
    histogram_new1.SetMaximum(2.0)
    histogram_new2.SetMaximum(2.0)
    
    
    #histogram_old.SetLineColor(2)
    histogram_new1.SetLineColor(3)
    histogram_new2.SetLineColor(4)

    #histogram_old.SetMarkerColor(2)
    histogram_new1.SetMarkerColor(3)
    histogram_new2.SetMarkerColor(4)

    #histogram_old.SetMarkerStyle(20)
    histogram_new1.SetMarkerStyle(21)
    histogram_new2.SetMarkerStyle(22)

    #histogram_old.SetMarkerSize(0.7)
    histogram_new1.SetMarkerSize(0.7)
    histogram_new2.SetMarkerSize(0.7)

    #histogram_old.SetFillColor(2)
    histogram_new1.SetFillColor(3)
    histogram_new2.SetFillColor(4)

    #histogram_old.SetFillStyle(3004)
    histogram_new1.SetFillStyle(3005)
    histogram_new2.SetFillStyle(3002)

    #xAxis_top = histogram_old.GetXaxis()
    xAxis_top = histogram_new1.GetXaxis()
    xAxis_top.SetTitle("Tau p_{T} (GeV)");
    xAxis_top.SetTitleOffset(1.2);
    xAxis_top.SetTitleSize(0.03)
    xAxis_top.SetLabelSize(0.03)


    #yAxis_top = histogram_old.GetYaxis()
    yAxis_top = histogram_new1.GetYaxis()
    yAxis_top.SetTitle("Data/MC Sf")
    yAxis_top.SetTitleOffset(1.2)
    yAxis_top.SetTitleSize(0.03)
    yAxis_top.SetLabelSize(0.03)
    yAxis_top.SetTickLength(0.04)
    
    legendTextSize = 0.040
    legendPosX     = 0.740
    legendPosY     = 0.510
    legendSizeX    = 0.190
    legendSizeY    = 0.420

    legend = ROOT.TLegend(0.7, 0.7, 0.85, 0.85, "", "brNDC")
    legend.SetFillStyle(0)
    legend.SetFillColor(10)
    legend.SetTextSize(0.018)
    #legend.AddEntry(histogram_old, "old", "f")
    legend.AddEntry(histogram_new1, "(Data - Bg)/ZTT", "f")
    legend.AddEntry(histogram_new2, "Data/(ZTT + Bg)", "f")

    #histogram_old.GetXaxis().SetRangeUser(20.,200.) ## Konstantin's old plots go till 1000 GeV but ours only till 200 GeV
    histogram_new1.GetXaxis().SetRangeUser(20.,200.)
    histogram_new2.GetXaxis().SetRangeUser(20.,200.)

    #histogram_old.GetYaxis().SetRangeUser(0., 1.4)
    histogram_new1.GetYaxis().SetRangeUser(0., 1.4)
    histogram_new2.GetYaxis().SetRangeUser(0., 1.4)

    ## --- PLOTTING ON THE SAME CANVAS ---##
    #histogram_old.Draw('E2') 
    histogram_new1.Draw('E2 same') 
    histogram_new2.Draw('E2 same') 
    
    legend.Draw()

    canvas.Update()
    canvas.Print(outputFileName + ".pdf")
    canvas.Print(outputFileName + ".png")
    canvas.Print(outputFileName + ".root")


for era in Eras:
    #FullInputFilePath_old = "{}/{}_tauTriggerEff_DeepTau2017v2p1.root".format(InputFilePath_old, era)
    FullInputFilePath_new_sub_from_data = '{}/turn_on_{}_subtract-from-data_fitted_LATEST.root'.format(InputFilePath_new, era)
    FullInputFilePath_new_add_to_dy_mc = '{}/turn_on_{}_add-to-dy-mc_fitted_LATEST.root'.format(InputFilePath_new, era)
    FullInputFilePath_new_TriggerSFs = '{}/NewTriggerSFs_{}.root'.format(InputFilePath_new, era)
    #f_old = ROOT.TFile.Open(FullInputFilePath_old, "READ")
    f_new_sub_from_data = ROOT.TFile.Open(FullInputFilePath_new_sub_from_data, "READ")
    f_new_add_to_dy_mc = ROOT.TFile.Open(FullInputFilePath_new_add_to_dy_mc, "READ")
    f_new_TriggerSFs = ROOT.TFile.Open(FullInputFilePath_new_TriggerSFs, "READ")
    for dm in Decay_modes:
        for chn in Channels:
            for wp in Working_points:
                histo_dict_OldSFs = {}
                histo_dict_NewSFs = {}
                histoName = "sf_{}_{}_dm{}_fitted".format(chn, wp, dm)
                histoName_sf_true = "{}_{}_dm{}_fitted_sf_true".format(chn, wp, dm)
                histoName_sf_fake = "{}_{}_dm{}_fitted_sf_fake".format(chn, wp, dm)
                #h_old = f_old.Get(histoName)
                h_new_sub_from_data = f_new_sub_from_data.Get(histoName)
                h_new_add_to_dy_mc = f_new_add_to_dy_mc.Get(histoName)
                h_sf_true = f_new_TriggerSFs.Get(histoName_sf_true)
                h_sf_fake = f_new_TriggerSFs.Get(histoName_sf_fake)
                outFileName_OldSFs = "{}/sf_{}_{}_{}_dm{}_fitted".format(OutputFilePath, era, chn, wp, dm)
                outFileName_NewSFs = "{}/New_sf_{}_{}_{}_dm{}_fitted".format(OutputFilePath, era, chn, wp, dm)
                #histo_dict_OldSFs["old"] = h_old
                histo_dict_OldSFs["sub_from_data"] = h_new_sub_from_data
                histo_dict_OldSFs["add_to_dy_mc"] = h_new_add_to_dy_mc
                histo_dict_NewSFs["sf_true"] = h_sf_true
                histo_dict_NewSFs["sf_fake"] = h_sf_fake
                makePlot_OldSFs(histo_dict_OldSFs, outFileName_OldSFs)
                makePlot_NewSFs(histo_dict_NewSFs, outFileName_NewSFs)
    #f_old.Close()
    f_new_sub_from_data.Close()
    f_new_add_to_dy_mc.Close()
