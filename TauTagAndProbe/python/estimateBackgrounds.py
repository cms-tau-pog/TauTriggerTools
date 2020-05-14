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

path_prefix = '' if 'TauTriggerTools' in os.getcwd() else 'TauTriggerTools/'
sys.path.insert(0, path_prefix + 'Common/python')
from AnalysisTools import *
#ROOT.ROOT.EnableImplicitMT(4)
ROOT.gROOT.SetBatch(True)
ROOT.gInterpreter.Declare('#include "{}TauTagAndProbe/interface/PyInterface.h"'.format(path_prefix))

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

# step 1: determine scale-factor for W+jets background in SS region
df_SS_high_mT = df_input.Filter("selection == %i" % selection_SS_high_mT)
sum_SS_high_mT = {}
for process in processes:
    df_SS_high_mT_process = df_SS_high_mT.Filter("type == %i" % get_type(process))
    sum_SS_high_mT[process] = df_SS_high_mT_process.Sum("weight")
    print("sum_SS_high_mT['%s'] = %1.2f (%i)" % (process, sum_SS_high_mT[process].GetValue(), df_SS_high_mT_process.Count().GetValue()))
sf_w_mc_SS = (sum_SS_high_mT['data'].GetValue() - (sum_SS_high_mT["ztt-mc"].GetValue() + sum_SS_high_mT['zmm-mc'].GetValue() + sum_SS_high_mT['ttbar-mc'].GetValue()))/sum_SS_high_mT['w-mc'].GetValue()
print("sf_w_mc_SS = %1.2f" % sf_w_mc_SS)

# step 2: determine QCD multijet background in SS region 
#        (Note: QCD multijet background in SS high mT sideband assumed to be negligible)
df_SS_low_mT = df_input.Filter("selection == %i" % selection_SS_low_mT)
sum_SS_low_mT = {}
for process in processes:
    df_SS_low_mT_process = df_SS_low_mT.Filter("type == %i" % get_type(process))
    sum_SS_low_mT[process] = df_SS_low_mT_process.Sum("weight")
    print("sum_SS_low_mT['%s'] = %1.2f (%i)" % (process, sum_SS_low_mT[process].GetValue(), df_SS_low_mT_process.Count().GetValue()))

# define SS->OS extrapolation factor for QCD multijet background 
sf_qcd_SS_to_OS = 1.

# step 3: determine scale-factor for W+jets background in OS region 
#        (Note: QCD multijet background in OS high mT sideband assumed to be negligible)
df_OS_high_mT = df_input.Filter("selection == %i" % selection_OS_high_mT)
sum_OS_high_mT = {}
for process in processes:
    df_OS_high_mT_process = df_OS_high_mT.Filter("type == %i" % get_type(process))
    sum_OS_high_mT[process] = df_OS_high_mT_process.Sum("weight")
    print("sum_OS_high_mT['%s'] = %1.2f (%i)" % (process, sum_OS_high_mT[process].GetValue(), df_OS_high_mT_process.Count().GetValue()))
sf_w_mc_OS = (sum_OS_high_mT['data'].GetValue() - (sum_OS_high_mT["ztt-mc"].GetValue() + sum_OS_high_mT['zmm-mc'].GetValue() + sum_OS_high_mT['ttbar-mc'].GetValue()))/sum_OS_high_mT['w-mc'].GetValue()
print("sf_w_mc_OS = %1.2f" % sf_w_mc_OS)

# step 4: print event yields in "signal" region for input RDataFrame objects
df_OS_low_mT = df_input.Filter("selection == %i" % selection_OS_low_mT)
sum_OS_low_mT = {}
for process in processes:
    df_OS_low_mT_process = df_OS_low_mT.Filter("type == %i" % get_type(process))
    sum_OS_low_mT[process] = df_OS_low_mT_process.Sum("weight")
    print("sum_OS_low_mT['%s'] = %1.2f (%i)" % (process, sum_OS_low_mT[process].GetValue(), df_OS_low_mT_process.Count().GetValue()))

# step 5: build RDataFrame object for 'data'
final_weight_data = ROOT.final_weight_data.Initialize(sf_qcd_SS_to_OS, sf_w_mc_OS, sf_w_mc_SS)
if args.mode == "subtract-from-data":
    df_output_data = df_input.Filter("(selection == %i && type == %i) || (selection == %i) || (selection == %i && (type == %i || type == %i || type == %i))" % (selection_OS_low_mT, type_data, selection_SS_low_mT, selection_OS_low_mT, type_zmm_mc, type_w_mc, type_ttbar_mc))
    df_output_data.Define('final_weight', "final_weight_data::GetDefault().operator()(selection, type, weight)")
elif args.mode == "add-to-dy-mc":
    df_output_data = df_input.Filter("selection == %i && type == %i" % (selection_OS_low_mT, type_data))

# step 6: build RDataFrame object for 'dy-mc'
final_weight_dy_mc = ROOT.final_weight_dy_mc.Initialize(sf_qcd_SS_to_OS, sf_w_mc_OS, sf_w_mc_SS)
if args.mode == "subtract-from-data":
    df_output_dy_mc = df_input.Filter("selection == %i && type == %i" % (selection_OS_low_mT, type_ztt_mc))
elif args.mode == "add-to-dy-mc":
    df_output_dy_mc = df_input.Filter("(selection == %i && type == %i) || (selection == %i) || (selection == %i && (type == %i || type == %i || type == %i)" % (selection_OS_low_mT, type_ztt_mc, selection_SS_low_mT, selection_OS_low_mT, type_zmm-mc, type_w_mc, type_ttbar_mc))
    df_output_dy_mc.Define('final_weight', "final_weight_dy_mc::GetDefault().operator()(selection, type, weight)")

# step 7: print data and dy-mc event yields in "signal" region for output RDataFrame objects
print "df_output_data = ", df_output_data
sum_OS_low_mT_data = df_output_data.Sum("weight")
print("sum_OS_low_mT_data = %1.2f (%i)" % (sum_OS_low_mT_data.GetValue(), df_output_data.Count().GetValue()))
print "df_output_dy_mc = ", df_output_dy_mc
sum_OS_low_mT_dy_mc = df_output_dy_mc.Sum("weight")
print("sum_OS_low_mT_dy_mc = %1.2f (%i)" % (sum_OS_low_mT_dy_mc.GetValue(), df_output_dy_mc.Count().GetValue()))

# step 8: write RDataFrame objects to output files
df_output_data.Snapshot('events', args.output_data)
df_output_dy_mc.Snapshot('events', args.output_dy_mc)
