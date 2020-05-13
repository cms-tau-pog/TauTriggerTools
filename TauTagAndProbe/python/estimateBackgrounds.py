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

if not(args.mode == "subtract-from-data" or args.mode == "add-to-dy-mc"):
    raise ValueError("Invalid configuration parameter mode = '%s' !!" % args.mode)

ROOT.gROOT.SetBatch(True)


path_prefix = '' if 'TauTriggerTools' in os.getcwd() else 'TauTriggerTools/'
sys.path.insert(0, path_prefix + 'Common/python')
from AnalysisTools import *


input_vec = ListToStdVector(args.input)

df_input = ROOT.RDataFrame('events', input_vec)

df_output_data  = None
df_output_dy_mc = None

processes = [ "data", "ztt-mc", "zmm-mc", "w-mc", "ttbar-mc" ]

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
#----------------------------------------------------------------------------------------------------

def get_type(process):
    if process == "data":
        return process_data
    elif process == "ztt-mc":
        process_ztt_mc
    elif process == "zmm-mc":
        process_zmm_mc
    elif process == "w-mc":
        process_w_mc
    elif process == "ttbar-mc":
        process_ttbar_mc
    else:
        raise ValueError("Invalid function argument: process = '%s' !!" % process)

# step 1: determine scale-factor for W+jets background in SS region
df_SS_high_mT = df_input.Filter("selection == %i" % selection_SS_high_mT)
sum_SS_high_mT = {}
for process in processes:
    sum_SS_high_mT[process] = df_SS_high_mT.Filter("type == %i" % get_type(process)).Sum("weight")
    print("sum_SS_high_mT['%s'] = %1.2f" % (process, sum_SS_high_mT[process]))
sf_w_mc_SS = (sum_SS_high_mT['data'] - (sum_SS_high_mT["ztt-mc"] + sum_SS_high_mT['zmm-mc'] + sum_SS_high_mT['ttbar-mc']))/sum_SS_high_mT['w-mc']
print("sf_w_mc_SS = %1.2f" % sf_w_mc_SS)

# step 2: determine QCD multijet background in SS region 
#        (Note: QCD multijet background in SS high mT sideband assumed to be negligible)
df_SS_low_mT = df_input.Filter("selection == %i" % selection_SS_low_mT)
for process in processes:
    sum_SS_low_mT[process] = df_SS_low_mT.Filter("type == %i" % get_type(process)).Sum("weight")
    print("sum_SS_low_mT['%s'] = %1.2f" % (process, sum_SS_low_mT[process]))

# define SS->OS extrapolation factor for QCD multijet background 
sf_qcd_SS_to_OS = 1.

# step 3: determine scale-factor for W+jets background in OS region 
#        (Note: QCD multijet background in OS high mT sideband assumed to be negligible)
df_OS_high_mT = df_input.Filter("selection == %i" % selection_OS_high_mT)
sum_OS_high_mT = {}
for process in processes:
    sum_OS_high_mT[process] = df_OS_high_mT.Filter("type == %i" % get_type(process)).Sum("weight")
    print("sum_OS_high_mT['%s'] = %1.2f" % (process, sum_OS_high_mT[process]))
sf_w_mc_OS = (sum_OS_high_mT['data'] - (sum_OS_high_mT["ztt-mc"] + sum_OS_high_mT['zmm-mc'] + sum_OS_high_mT['ttbar-mc']))/sum_OS_high_mT['w-mc']
print("sf_w_mc_OS = %1.2f" % sf_w_mc_OS)

# step 4: build RDataFrame object for 'data'
def final_weight_data(selection, type, weight):
    final_weight = None
    if selection == selection_OS_low_mT and type == type_data:
        final_weight = 1.
    elif selection == selection_SS_low_mT:
        final_weight = sf_qcd_SS_to_OS
        if type == 'data':
            final_weight *= -1. 
        else:
            final_weight *= +1. * weight
            if type == 'w-mc':
                final_weight *= sf_w_mc_SS
    elif selection == selection_OS_low_mT and (type == type_zmm_mc or type == type_w_mc or type == type_ttbar_mc):
        final_weight = -1. * weight
        if type == 'w-mc':
                final_weight *= sf_w_mc_OS
    else:
        raise ValueError("Invalid function arguments: selection = '%i', type = '%i' !!" % (selection, type))

if args.mode == "subtract-from-data":
    df_output_data = df_input.Filter("(selection == %i && type == %i) || (selection == %i) || (selection == %i && (type == %i || type == %i || type == %i))" % (selection_OS_low_mT, type_data, selection_SS_low_mT, selection_OS_low_mT, type_zmm_mc, type_w_mc, type_ttbar_mc))
    df_output_data.Define('final_weight', final_weight_data, { "selection", "type", "weight" })
elif args.mode == "add-to-dy-mc":
    df_output_data = df_input.Filter("selection == %i && type == %i" % (selection_OS_low_mT, type_data))

# step 5: build RDataFrame object for 'dy-mc'
def final_weight_dy_mc(selection, type, weight):
    final_weight = None
    if selection == selection_OS_low_mT and type == type_ztt_mc:
        final_weight = +1. * weight
    elif selection == selection_SS_low_mT:
        final_weight = sf_qcd_SS_to_OS
        if type == 'data':
            final_weight *= +1. 
        else:
            final_weight *= -1. * weight
            if type == type_w_mc:
                final_weight *= sf_w_mc_SS
    elif selection == selection_OS_low_mT and (type == type_zmm_mc or type == type_w_mc or type == type_ttbar_mc):
        final_weight = +1. * weight
        if type == 'w-mc':
                final_weight *= sf_w_mc_OS
    else:
        raise ValueError("Invalid function arguments: selection = '%i', type = '%i' !!" % (selection, type))

if args.mode == "subtract-from-data":
    df_output_dy_mc = df_input.Filter("selection == %i && type == %i" % (selection_OS_low_mT, type_ztt_mc))
elif args.mode == "add-to-dy-mc":
    df_output_dy_mc = df_input.Filter("(selection == %i && type == %i) || (selection == %i) || (selection == %i && (type == %i || type == %i || type == %i)" % (selection_OS_low_mT, type_ztt_mc, selection_SS_low_mT, selection_OS_low_mT, type_zmm-mc, type_w_mc, type_ttbar_mc))
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
