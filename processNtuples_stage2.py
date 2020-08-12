import os, subprocess, sys

Area = "/hdfs/local/ram/run2_ntuples" 

Eras = ["2016", "2017", "2018"]
Sidebands = ["signal", "w_enriched", "OS_low_mT", "SS_low_mT", "OS_high_mT", "SS_high_mT"]
Types = ["data", "ztt_mc", "zmm_mc", "w_mc", "ttbar_mc"]
Input_File_Label = ["SingleMuon", "DYJetsToLL", "WJetsToLNu", "TTTo2L2Nu", "TTToSemiLeptonic"]
Output_File_Label = ["DATA", "DY", "W", "TT_2L", "TT_SL"]

Lumiscale = { ## Defined as (x-sec * Integ.Lumi)/Nevt
    "2016": {
        "DATA": 1.0,
        "DY": 1.49,
        "W": 25.41,
        "TT_2L": 0.047,
        "TT_SL": 0.122
    },
    "2017": {
        "DATA": 1.0,
        "DY": 2.58,
        "W": 56.89,
        "TT_2L": 0.053,
        "TT_SL": 0.13
    },
    "2018": {
        "DATA": 1.0,
        "DY": 3.66,
        "W": 51.71,
        "TT_2L": 0.082,
        "TT_SL": 0.214
    }
}



Lumiscale_num = { ## Defined as (x-sec * Integ.Lumi)
    "2016": {
        "DATA": 1.0,
        "DY": 218172198.0,
        "W": 2208808530.0,
        "TT_2L": 3173560.0,
        "TT_SL": 13122168.0 
    },
    "2017": {
        "DATA": 1.0,
        "DY": 252204630.0,
        "W": 2553358050.0,
        "TT_2L": 3668600.0,
        "TT_SL": 15169080.0 
    },
    "2018": {
        "DATA": 1.0,
        "DY": 362810034.0,
        "W": 3673143990.0,
        "TT_2L": 5277480.0,
        "TT_SL": 21821544.0 
    }
}

N_total = { ## Total number of unweighted MC events
    "2016": {
        "DATA": 1.0,
        "DY": 146280395.0,
        "W": 86916455.0, 
        "TT_2L": 67926800.0,        
        "TT_SL": 107604800.0 
    },
    "2017": {
        "DATA": 1.0,
        "DY": 97800939.0, 
        "W": 44881137.0,
        "TT_2L": 69155808.0,
        "TT_SL": 110014744.0 
    },
    "2018": {
        "DATA": 1.0,
        "DY": 99100315.0,
        "W": 71026861.0,
        "TT_2L": 64310000.0,
        "TT_SL": 101550000.0 
    }
}





def run_cmd(command):
  print "executing command = '%s'" % command
  p = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  stdout, stderr = p.communicate()
  return stdout

print("Creating directories to store the skimmed Ntuples (if they are not already created)")
run_cmd('mkdir -p tuples')
run_cmd('mkdir -p tuples/skimmed-2016')
run_cmd('mkdir -p tuples/skimmed-2017')
run_cmd('mkdir -p tuples/skimmed-2018')
print("Clearing up old stuff (if present)")
run_cmd('rm tuples/skimmed-201*/*')


for era in Eras:
  for sample in Types:
    for sideband in Sidebands:
      print("Era: {}, Type: {}, Sideband: {}".format(era, sample, sideband))
      if(sample == "data"):
        ID = 0
        print("Running skimTuple_NEW.py for ID: {}, Input_File_Label[ID]: {}, Output_File_Label[ID]: {}".format(str(ID), Input_File_Label[ID], Output_File_Label[ID]))
        run_cmd("python TauTagAndProbe/python/skimTuple.py --input {}/full/{}/{}* --config TauTagAndProbe/data/{}/triggers.json --selection {} --output tuples/skimmed-{}/{}-{}-{}.root --type {} --sideband {} --lumiScale {}".format( Area, era, Input_File_Label[ID], era, "DeepTau", era, Output_File_Label[ID], sample, sideband, sample, sideband, str(Lumiscale[era][Output_File_Label[ID]]) ))
      elif((sample == "ztt_mc") or (sample == "zmm_mc")):
        ID = 1
        print("Running skimTuple_NEW.py for ID: {}, Input_File_Label[ID]: {}, Output_File_Label[ID]: {}".format(str(ID), Input_File_Label[ID], Output_File_Label[ID]))
        run_cmd("python TauTagAndProbe/python/skimTuple.py --input {}/full/{}/{}* --config TauTagAndProbe/data/{}/triggers.json --selection {} --output tuples/skimmed-{}/{}-{}-{}.root --type {} --pu {}/Pileup_Data{}.root --sideband {} --lumiScale {}".format( Area, era, Input_File_Label[ID], era, "DeepTau", era, Output_File_Label[ID], sample, sideband, sample, Area, era, sideband, str(Lumiscale[era][Output_File_Label[ID]]) ))
      elif(sample == "w_mc"):
        ID = 2
        print("Running skimTuple_NEW.py for ID: {}, Input_File_Label[ID]: {}, Output_File_Label[ID]: {}".format(str(ID), Input_File_Label[ID], Output_File_Label[ID]))
        run_cmd("python TauTagAndProbe/python/skimTuple.py --input {}/full/{}/{}* --config TauTagAndProbe/data/{}/triggers.json --selection {} --output tuples/skimmed-{}/{}-{}-{}.root --type {} --pu {}/Pileup_Data{}.root --sideband {} --lumiScale {}".format( Area, era, Input_File_Label[ID], era, "DeepTau", era, Output_File_Label[ID], sample, sideband, sample, Area, era, sideband, str(Lumiscale[era][Output_File_Label[ID]]) ))
      elif(sample == "ttbar_mc"):
        ID = 3  
        print("Running skimTuple_NEW.py for ID: {}, Input_File_Label[ID]: {}, Output_File_Label[ID]: {}".format(str(ID), Input_File_Label[ID], Output_File_Label[ID]))
        run_cmd("python TauTagAndProbe/python/skimTuple.py --input {}/full/{}/{}* --config TauTagAndProbe/data/{}/triggers.json --selection {} --output tuples/skimmed-{}/{}-{}-{}.root --type {} --pu {}/Pileup_Data{}.root --sideband {} --lumiScale {}".format( Area, era, Input_File_Label[ID], era, "DeepTau", era, Output_File_Label[ID], sample, sideband, sample, Area, era, sideband, str(Lumiscale[era][Output_File_Label[ID]]) ))
        ID = 4  
        print("Running skimTuple_NEW.py for ID: {}, Input_File_Label[ID]: {}, Output_File_Label[ID]: {}".format(str(ID), Input_File_Label[ID], Output_File_Label[ID]))
        run_cmd("python TauTagAndProbe/python/skimTuple.py --input {}/full/{}/{}* --config TauTagAndProbe/data/{}/triggers.json --selection {} --output tuples/skimmed-{}/{}-{}-{}.root --type {} --pu {}/Pileup_Data{}.root --sideband {} --lumiScale {}".format( Area, era, Input_File_Label[ID], era, "DeepTau", era, Output_File_Label[ID], sample, sideband, sample, Area, era, sideband, str(Lumiscale[era][Output_File_Label[ID]]) ))

  print("Adding all root files for era {}".format(era))        
  run_cmd("hadd -f input_stage2p5_{}.root tuples/skimmed-{}/*.root".format(era, era))
  

  print("Creating input files for CreateTrunOn.py script for era {} and mode {}".format(era, "subtract-from-data"))        
  run_cmd("python TauTagAndProbe/python/estimateBackgrounds.py --input input_stage2p5_{}.root --mode {} --output-data estimateBackgrounds_{}_DATA_mode_{} --output-dy-mc estimateBackgrounds_{}_DY_MC_mode_{} --output-w-enriched w_enriched_{} --output-signal signal_{}".format(era, "subtract-from-data", era, "subtract-from-data",  era, "subtract-from-data", era, era))

  print("Creating input files for CreateTrunOn.py script for era {} and mode {}".format(era, "add-to-dy-mc"))        
  run_cmd("python TauTagAndProbe/python/estimateBackgrounds.py --input input_stage2p5_{}.root --mode {} --output-data estimateBackgrounds_{}_DATA_mode_{} --output-dy-mc estimateBackgrounds_{}_DY_MC_mode_{}".format(era, "add-to-dy-mc", era, "add-to-dy-mc",  era, "add-to-dy-mc"))
  
  print("Running CreateTrunOn.py script for era {} and {} region".format(era, "signal"))
  run_cmd("python TauTagAndProbe/python/createTrunOn.py --input-data signal_{}_data.root --input-dy-mc signal_{}_ztt_mc.root --output turn_on_{}_signal --branchname-weight-data {} --branchname-weight-dy-mc {}".format(era, era, era, "weight", "weight"))

  print("Running CreateTrunOn.py script for era {} and {} region".format(era, "w-enriched"))
  run_cmd("python TauTagAndProbe/python/createTrunOn.py --input-data  w_enriched_{}_data.root --input-dy-mc  w_enriched_{}_mc.root --output turn_on_{}_w_enriched --branchname-weight-data {} --branchname-weight-dy-mc {}".format(era, era, era, "weight", "weight"))

  print("Running CreateTrunOn.py script for era {} and mode {}".format(era, "subtract-from-data"))        
  run_cmd("python TauTagAndProbe/python/createTrunOn.py --input-data estimateBackgrounds_{}_DATA_mode_{}.root --input-dy-mc estimateBackgrounds_{}_DY_MC_mode_{}.root --output turn_on_{}_{}_LATEST --branchname-weight-data {} --branchname-weight-dy-mc {} ".format(era, "subtract-from-data", era, "subtract-from-data", era, "subtract-from-data", "final_weight", "weight"))
  
  print("Running CreateTrunOn.py script for era {} and mode {}".format(era, "add-to-dy-mc"))        
  run_cmd("python TauTagAndProbe/python/createTrunOn.py --input-data estimateBackgrounds_{}_DATA_mode_{}.root --input-dy-mc estimateBackgrounds_{}_DY_MC_mode_{}.root --output turn_on_{}_{}_LATEST --branchname-weight-data {} --branchname-weight-dy-mc {} ".format(era, "add-to-dy-mc", era, "add-to-dy-mc", era, "add-to-dy-mc", "weight", "final_weight"))

  print("Running fitTurnOn.py script for era {} and {} region".format(era, "signal"))        
  run_cmd("python TauTagAndProbe/python/fitTurnOn.py --input turn_on_{}_{}.root --output turn_on_{}_{}_fitted --mode {}".format(era, "signal", era, "signal", "fixed"))

  print("Running fitTurnOn.py script for era {} and {} region".format(era, "w_enriched"))        
  run_cmd("python TauTagAndProbe/python/fitTurnOn.py --input turn_on_{}_{}.root --output turn_on_{}_{}_fitted --mode {}".format(era, "w_enriched", era, "w_enriched", "fixed"))

  print("Running fitTurnOn.py script for era {} and mode {}".format(era, "subtract-from-data"))        
  run_cmd("python TauTagAndProbe/python/fitTurnOn.py --input turn_on_{}_{}_LATEST.root --output turn_on_{}_{}_fitted_LATEST --mode {}".format(era, "subtract-from-data", era, "subtract-from-data", "adaptive"))

  print("Running fitTurnOn.py script for era {} and mode {}".format(era, "add-to-dy-mc"))        
  run_cmd("python TauTagAndProbe/python/fitTurnOn.py --input turn_on_{}_{}_LATEST.root --output turn_on_{}_{}_fitted_LATEST --mode {}".format(era, "add-to-dy-mc", era, "add-to-dy-mc", "adaptive"))
  
  print("Running computeTriggerSFs.py script for era {}".format(era))
  run_cmd("python TauTagAndProbe/python/computeTriggerSFs.py --input_stage2p5 input_stage2p5_{}.root --input_signal turn_on_{}_signal_fitted.root --input_w_enriched turn_on_{}_w_enriched_fitted.root --output NewTriggerSFs_{}".format(era, era, era, era))

  print("Making directory {} to store the final plots and SF root files".format("Tau_Trigger_sf_plots"))
  run_cmd("mkdir -p $PWD/{}".format("Tau_Trigger_sf_plots"))
  print("Running TriggerSF_plotter.py script for era {}".format(era))
  run_cmd("python TauTagAndProbe/python/TriggerSF_plotter.py --era {} --inputFilePath-new $PWD --outputFilePath $PWD/{}".format(era, "Tau_Trigger_sf_plots"))
