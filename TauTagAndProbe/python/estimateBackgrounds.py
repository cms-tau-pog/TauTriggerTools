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
parser.add_argument('--mode', required=True, type=str, help="subtract backgrounds from data or add backgrounds to ZTT MC")
args = parser.parse_args()

if not(mode == "subtract-from-data" or mode == "add-to-dy-mc"):
    raise ValueError("Invalid configuration parameter mode = '%s' !!" % mode)

df_input = ROOT.RDataFrame('events', args.input)

df_output_data = ROOT.RDataFrame('events')
df_output_dy_mc = ROOT.RDataFrame('events')

processes = [ "data", "ztt-mc", "zmm-mc", "w-mc", "ttbar-mc" ]

# step 1: determine scale-factor for W+jets background in SS region
df_SS_high_mT = df_input.Filter("selection == SS_high_mT")
sum_SS_high_mT = {}
for process in processes:
    sum_SS_high_mT[process] = df_SS_high_mT.Filter("type == %s" % process).Sum("weight")
    print("sum_SS_high_mT['%s'] = %1.2f" % (process, sum_SS_high_mT[process]))
sf_w_mc_SS = (sum_SS_high_mT['data'] - (sum_SS_high_mT["ztt-mc"] + sum_SS_high_mT['zmm-mc'] + sum_SS_high_mT['ttbar-mc']))/sum_SS_high_mT['w-mc']
print("sf_w_mc_SS = %1.2f" % sf_w_mc_SS)

# step 2: determine QCD multijet background in SS region 
#        (Note: QCD multijet background in SS high mT sideband assumed to be negligible)
df_SS_low_mT = df_input.Filter("selection == SS_low_mT")
for process in processes:
    sum_SS_low_mT[process] = df_SS_low_mT.Filter("type == %s" % process).Sum("weight")
    print("sum_SS_low_mT['%s'] = %1.2f" % (process, sum_SS_low_mT[process]))

# define SS->OS extrapolation factor for QCD multijet background 
sf_qcd_SS_to_OS = 1.

# step 3: determine scale-factor for W+jets background in OS region 
#        (Note: QCD multijet background in OS high mT sideband assumed to be negligible)
df_OS_high_mT = df_input.Filter("selection == OS_high_mT")
sum_OS_high_mT = {}
for process in processes:
    sum_OS_high_mT[process] = df_OS_high_mT.Filter("type == %s" % process).Sum("weight")
    print("sum_OS_high_mT['%s'] = %1.2f" % (process, sum_OS_high_mT[process]))
sf_w_mc_OS = (sum_OS_high_mT['data'] - (sum_OS_high_mT["ztt-mc"] + sum_OS_high_mT['zmm-mc'] + sum_OS_high_mT['ttbar-mc']))/sum_OS_high_mT['w-mc']
print("sf_w_mc_OS = %1.2f" % sf_w_mc_OS)

# step 4: build RDataFrame object for 'data'
def final_weight_data(selection, type, weight):
    final_weight = None
    if selection == 'OS_low_mT' and type == 'data':
        final_weight = 1.
    elif selection == 'SS_low_mT':
        final_weight = sf_qcd_SS_to_OS
        if type == 'data':
            final_weight *= -1. 
        else:
            final_weight *= +1. * weight
            if type == 'w-mc':
                final_weight *= sf_w_mc_SS
    elif selection == 'OS_low_mT' and (type == 'zmm-mc' or type == 'w-jets-mc' or type == 'ttbar-mc'):
        final_weight = -1. * weight
        if type == 'w-mc':
                final_weight *= sf_w_mc_OS
    else:
        raise ValueError("Invalid function arguments: selection = '%s', type = '%s' !!" % (selection, type))

if mode == "subtract-from-data":
    df_output_data = df_input.Filter("(selection == 'OS_low_mT' && type == 'data') || (selection == 'SS_low_mT') || (selection == 'OS_low_mT' && (type == 'zmm-mc' || type == 'w-jets-mc' || type == 'ttbar-mc'))")
    df_output_data.Define('final_weight', final_weight_data, { "selection", "type", "weight" })
elif mode == "add-to-dy-mc":
    df_output_data = df_input.Filter("selection == 'OS_low_mT' && type == 'data'")

# step 5: build RDataFrame object for 'dy-mc'
def final_weight_dy_mc(selection, type, weight):
    final_weight = None
    if selection == 'OS_low_mT' and type == 'ztt-mc':
        final_weight = +1. * weight
    elif selection == 'SS_low_mT':
        final_weight = sf_qcd_SS_to_OS
        if type == 'data':
            final_weight *= +1. 
        else:
            final_weight *= -1. * weight
            if type == 'w-mc':
                final_weight *= sf_w_mc_SS
    elif selection == 'OS_low_mT' and (type == 'zmm-mc' or type == 'w-jets-mc' or type == 'ttbar-mc'):
        final_weight = +1. * weight
        if type == 'w-mc':
                final_weight *= sf_w_mc_OS
    else:
        raise ValueError("Invalid function arguments: selection = '%s', type = '%s' !!" % (selection, type))

if mode == "subtract-from-data":
    df_output_dy_mc = df_input.Filter("selection == 'OS_low_mT' && type == 'ztt-mc'")
elif mode == "add-to-dy-mc":
    df_output_dy_mc = df_input.Filter("(selection == 'OS_low_mT' && type == 'ztt-mc') || (selection == 'SS_low_mT') || (selection == 'OS_low_mT' && (type == 'zmm-mc' || type == 'w-jets-mc' || type == 'ttbar-mc')")
    df_output_dy_mc.Define('final_weight', final_weight_dy_mc, { "selection", "type", "weight" })

# step 6: write RDataFrame objects to output files
output_file_data = ROOT.TFile(args.output_data + '.root', 'RECREATE')
output_file_data.cd()
df_output_data.Write()
output_file_data.Close()

output_file_dy_mc = ROOT.TFile(args.output_dy_mc + '.root', 'RECREATE')
output_file_dy_mc.cd()
df_output_dy_mc.Write()
output_file_dy_mc.Close()
