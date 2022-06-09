#!/usr/bin/env python
import argparse
from array import array
import math
import numpy as np
import os
import re
import sys
import ROOT

parser = argparse.ArgumentParser(description='Skim full tuple.')
parser.add_argument('--input', required=True, type=str, nargs='+', help="input files")
parser.add_argument('--config', required=True, type=str, help="config with triggers description")
parser.add_argument('--selection', required=True, type=str, help="tau selection")
parser.add_argument('--output', required=True, type=str, help="output file")
parser.add_argument('--type', required=True, type=str, default='data', help="Define the sample type among the following: data, ztt_mc, zmm_mc, w_mc, ttbar_mc")
parser.add_argument('--lumiScale', required=True, type=float, default=1.0, help="LumiScale factor")
parser.add_argument('--sideband', required=True, type=str, default='signal', 
                    help="Event level selections to define sidebands: signal, w_enriched, OS_low_mT, OS_high_mT, SS_low_mT, SS_high_mT")
parser.add_argument('--pu', required=False, type=str, default=None,
                    help="file with the pileup profile for the data taking period")
args = parser.parse_args()

path_prefix = '' if 'TauTriggerTools' in os.getcwd() else 'TauTriggerTools/'
sys.path.insert(0, path_prefix + 'Common/python')
from AnalysisTypes import *
from AnalysisTools import *
import TriggerConfig
ROOT.ROOT.EnableImplicitMT(4)
ROOT.gROOT.SetBatch(True)
ROOT.gInterpreter.Declare('#include "{}TauTagAndProbe/interface/PyInterface.h"'.format(path_prefix))

if args.type not in ["data", "ztt_mc", "zmm_mc", "w_mc", "ttbar_mc"]:
    raise RuntimeError("Invalid sample type")

if args.sideband not in ['signal', 'w_enriched', 'OS_low_mT','OS_high_mT','SS_low_mT','SS_high_mT']:
    raise RuntimeError("Invalid sideband")
sideband = args.sideband
LumiScale = args.lumiScale
input_vec = ListToStdVector(args.input)

if args.type != 'data':
    if args.pu is None:
        raise RuntimeError("Pileup file should be provided for mc.")
    data_pu_file = ROOT.TFile(args.pu, 'READ')
    data_pu = data_pu_file.Get('pileup')
    df_all = ROOT.RDataFrame('all_events', input_vec)
    mc_pu = df_all.Histo1D(ROOT.RDF.TH1DModel(data_pu), 'npu')
    ROOT.PileUpWeightProvider.Initialize(data_pu, mc_pu.GetPtr())

trig_descriptors, channel_triggers = TriggerConfig.Load(args.config)
trigger_dict, filter_dict = TriggerConfig.LoadTriggerDictionary(input_vec)
triggerMatch = ROOT.TriggerMatchProvider.Initialize()
channels = {}
for channel_name, channel_trig_descs in channel_triggers.items():
    channel_id = ParseEnum(Channel, channel_name)
    channels[channel_name] = channel_id
    for desc in channel_trig_descs:
        if 'sample_types' in desc and args.type not in desc['sample_types']: continue
        if desc['leg_types'][-1] != 'tau': continue
        match_desc = ROOT.TriggerMatchProvider.MatchDescriptor()
        pattern = '^{}.*'.format(desc['name'])
        hlt_paths = TriggerConfig.GetMatchedTriggers(trigger_dict, pattern)
        match_desc.match_mask = int(TriggerConfig.GetMatchMask(hlt_paths))
        filter_names = desc['filters'][-1]
        match_desc.filter_hashes = ListToStdVector([ filter_dict[f] for f in filter_names ], elem_type='UInt_t')
        if 'min_run' in desc and args.type == 'data':
            match_desc.min_run = desc['min_run']
        if 'max_run' in desc and args.type == 'data':
            match_desc.max_run = desc['max_run']
        sel_name = 'selection_' + channel_name
        if sel_name in desc:
            if 'hltObj_pt' in desc[sel_name]:
                match_desc.hltObj_pt = desc[sel_name]['hltObj_pt']
            if 'l1Tau_pt' in desc[sel_name]:
                match_desc.l1Tau_pt = desc[sel_name]['l1Tau_pt']
            if 'l1Tau_hwIso' in desc[sel_name]:
                match_desc.l1Tau_hwIso = desc[sel_name]['l1Tau_hwIso']
        triggerMatch.Add(channel_id, match_desc)

selection_id = ParseEnum(TauSelection, args.selection)
df = ROOT.RDataFrame('events', input_vec)
process_id = ParseEnum(Process, args.type)
sideband_id = ParseEnum(SideBand, args.sideband)
df = df.Define('type', str(process_id))
df = df.Define('selection', str(sideband_id))

if((sideband == "signal") or (sideband == "OS_low_mT")):
    df = df.Filter('''
                   (tau_sel & {}) != 0 && muon_pt > 27 && muon_iso < 0.1 && muon_mt < 30
                   && tau_pt > 20 && abs(tau_eta) < 2.1 && tau_decayMode != 5 && tau_decayMode != 6
                   && vis_mass > 40 && vis_mass < 80 && muon_charge != tau_charge
                   '''.format(selection_id)) ## SIGNAL REGION
elif sideband == "w_enriched":
    df = df.Filter(''' 
                   (tau_sel & {}) != 0 && muon_pt > 27 && muon_iso < 0.1 && muon_mt > 50
                   && tau_pt > 20 && abs(tau_eta) < 2.1 && tau_decayMode != 5 && tau_decayMode != 6
                   '''.format(selection_id)) ## (W-ENRICHED SIDEBAND)
elif sideband == "SS_low_mT":
      df = df.Filter('''
                     (tau_sel & {}) != 0 && muon_pt > 27 && muon_iso < 0.1 && muon_mt < 30
                     && tau_pt > 20 && abs(tau_eta) < 2.1 && tau_decayMode != 5 && tau_decayMode != 6
                     && vis_mass > 40 && vis_mass < 80 && muon_charge == tau_charge
                     '''.format(selection_id))
elif sideband == "OS_high_mT":
      df = df.Filter('''
                      (tau_sel & {}) != 0 && muon_pt > 27 && muon_iso < 0.1 && muon_mt > 70 && muon_mt < 120
                      && tau_pt > 20 && abs(tau_eta) < 2.1 && tau_decayMode != 5 && tau_decayMode != 6
                      && vis_mass > 40 && vis_mass < 80 && muon_charge != tau_charge
                      '''.format(selection_id))
elif sideband == "SS_high_mT":
      df = df.Filter('''
                     (tau_sel & {}) != 0 && muon_pt > 27 && muon_iso < 0.1 && muon_mt > 70 && muon_mt < 120
                     && tau_pt > 20 && abs(tau_eta) < 2.1 && tau_decayMode != 5 && tau_decayMode != 6
                     && vis_mass > 40 && vis_mass < 80 && muon_charge == tau_charge
                     '''.format(selection_id))


if selection_id == TauSelection.DeepTau:
    df = df.Filter('(byDeepTau2017v2p1VSmu & (1 << {})) != 0'.format(DiscriminatorWP.Tight))

if args.type == 'data':
     df = df.Define('puWeight', "float(1.)")
     df = df.Define("genEventWeight_signOnly", "0.")
     df = df.Define('lumiScale', str(LumiScale)) 
     if sideband == "w_enriched":
         df = df.Define('weight', "muon_charge != tau_charge ? 1. : -1.") ## W-ENRICHED SIDEBAND FOR DATA
         print("w_enriched data yield ", df.Sum("weight").GetValue())
     else:
         df = df.Define('weight', "1.")
else:
     if sideband == "w_enriched":
         df = df.Filter('tau_charge + muon_charge == 0') ## W-ENRICHED SIDEBAND FOR MC
     if args.type == 'ztt_mc':
         print("Applying tau_gen_match == 5 cut to process ztt_mc in sideband {}".format(sideband))
         df = df.Filter('tau_gen_match == 5')
     elif( (args.type == 'zmm_mc') or 
           ( ((sideband == "w_enriched") or (sideband == "signal")) 
             and ((args.type == 'w_mc') or (args.type == 'ttbar_mc')) ) ):
         print("Applying tau_gen_match != 5 cut to process {} in sideband {}".format(args.type, sideband))
         df = df.Filter('tau_gen_match != 5') ## THIS CUT WASN'T APPLIED TO ANY OTHER "NON Z-MM" PROCESS IN ORIGINAL STUDY
     df = df.Define('puWeight', "PileUpWeightProvider::GetDefault().GetWeight(npu)")
     df = df.Define('genEventWeight_signOnly', "genEventWeight >= 0. ? +1. : -1.")
     N_eff = float(df.Sum("genEventWeight_signOnly").GetValue())
     N_tot = float(df.Count().GetValue()) 
     print("N_tot = %1.2f, N_eff = %1.2f" % (N_tot, N_eff))
     df = df.Define('lumiScale', "%s * %1.2f / %1.2f" % (str(LumiScale), N_tot, N_eff)) # LumiScale = x-sec * Integ. Lumi. 
     print("lumiScale = %1.2f" % eval("%s * %1.2f / %1.2f" % (str(LumiScale), N_tot, N_eff)))
     df = df.Define('weight', "puWeight * genEventWeight_signOnly * lumiScale")
 
skimmed_branches = [
    'type', 'selection', 
    'puWeight', 'genEventWeight', 'genEventWeight_signOnly', 'lumiScale', 'weight', 
    'tau_pt', 'tau_eta', 'tau_phi', 'tau_mass', 'tau_charge', 'tau_decayMode', 
    'byIsolationMVArun2017v2DBoldDMwLT2017', 'byDeepTau2017v2p1VSjet'
]

deltaRThr = 0.5
for channel_name, channel_id in channels.items():
    pass_branch = str('pass_' + channel_name)
    df = df.Define(pass_branch, '''TriggerMatchProvider::GetDefault().Pass({}, run, tau_eta, tau_phi, hlt_accept, {},
                   hltObj_types, hltObj_pt, hltObj_eta, hltObj_phi, hltObj_hasPathName, filter_hltObj, filter_hash,
                   l1Tau_pt, l1Tau_hwIso)'''.format(channel_id, deltaRThr))
    skimmed_branches.append(pass_branch)

df.Snapshot('events', args.output, ListToStdVector(skimmed_branches))
