##--- STEPS TO COMPUTE TRIGGER SCALE FACTORS ------##

(1) Setup the area following this recipe:

cmsrel CMSSW_10_2_20

cd CMSSW_10_2_20/src

cmsenv

git cms-init

git cms-addpkg RecoMET/METFilters

git cms-merge-topic cms-egamma:EgammaPostRecoTools

git clone https://github.com/kandrosov/TauTriggerTools $CMSSW_BASE/src/TauTriggerTools

cd $CMSSW_BASE/src/TauTriggerTools

git branch -a

git checkout remotes/origin/new-tuple

cd $CMSSW_BASE/src

scram b -j 10

(2) Submit CRAB jobs (after setting up CRAB credentials) for all the samples and all the eras via. the commands (change --site option to the site where you have write access):

    Common/scripts/crab_submit.py  --workArea work-area-2016 --cfg TauTagAndProbe/test/produceTuples.py --site T2_EE_Estonia --output tau_hlt_prod_2016 TauTagAndProbe/data/2016/crab/TT.txt
    Common/scripts/crab_submit.py  --workArea work-area-2017 --cfg TauTagAndProbe/test/produceTuples.py --site T2_EE_Estonia --output tau_hlt_prod_2017 TauTagAndProbe/data/2017/crab/TT.txt
    Common/scripts/crab_submit.py  --workArea work-area-2018 --cfg TauTagAndProbe/test/produceTuples.py --site T2_EE_Estonia --output tau_hlt_prod_2018 TauTagAndProbe/data/2018/crab/TT.txt
    Common/scripts/crab_submit.py  --workArea work-area-2016 --cfg TauTagAndProbe/test/produceTuples.py --site T2_EE_Estonia --output tau_hlt_prod_2016 TauTagAndProbe/data/2016/crab/W.txt
    Common/scripts/crab_submit.py  --workArea work-area-2017 --cfg TauTagAndProbe/test/produceTuples.py --site T2_EE_Estonia --output tau_hlt_prod_2017 TauTagAndProbe/data/2017/crab/W.txt
    Common/scripts/crab_submit.py  --workArea work-area-2018 --cfg TauTagAndProbe/test/produceTuples.py --site T2_EE_Estonia --output tau_hlt_prod_2018 TauTagAndProbe/data/2018/crab/W.txt
    Common/scripts/crab_submit.py  --workArea work-area-2016 --cfg TauTagAndProbe/test/produceTuples.py --site T2_EE_Estonia --output tau_hlt_prod_2016 TauTagAndProbe/data/2016/crab/DY.txt
    Common/scripts/crab_submit.py  --workArea work-area-2017 --cfg TauTagAndProbe/test/produceTuples.py --site T2_EE_Estonia --output tau_hlt_prod_2017 TauTagAndProbe/data/2017/crab/DY.txt
    Common/scripts/crab_submit.py  --workArea work-area-2018 --cfg TauTagAndProbe/test/produceTuples.py --site T2_EE_Estonia --output tau_hlt_prod_2018 TauTagAndProbe/data/2018/crab/DY.txt
    Common/scripts/crab_submit.py  --workArea work-area-2016 --cfg TauTagAndProbe/test/produceTuples.py --site T2_EE_Estonia --output tau_hlt_prod_2016 TauTagAndProbe/data/2016/crab/Data_SingleMuon.txt
    Common/scripts/crab_submit.py  --workArea work-area-2017 --cfg TauTagAndProbe/test/produceTuples.py --site T2_EE_Estonia --output tau_hlt_prod_2017 TauTagAndProbe/data/2017/crab/Data_SingleMuon.txt
    Common/scripts/crab_submit.py  --workArea work-area-2018 --cfg TauTagAndProbe/test/produceTuples.py --site T2_EE_Estonia --output tau_hlt_prod_2018 TauTagAndProbe/data/2018/crab/Data_SingleMuon_ABC.txt
    Common/scripts/crab_submit.py  --workArea work-area-2018 --cfg TauTagAndProbe/test/produceTuples.py --site T2_EE_Estonia --output tau_hlt_prod_2018 TauTagAndProbe/data/2018/crab/Data_SingleMuon_D.txt


(3) After job completion, merge the CRAB output files (for the respective samples and eras) into single root files using hadd and move them into a dedictated folder in your personal mass storage area (called "/hdfs/local/ram/run2_ntuples" below):


    mkdir -p /hdfs/local/ram/run2_ntuples
    mkdir -p /hdfs/local/ram/run2_ntuples/full
    mkdir -p /hdfs/local/ram/run2_ntuples/full/2016
    mkdir -p /hdfs/local/ram/run2_ntuples/full/2017
    mkdir -p /hdfs/local/ram/run2_ntuples/full/2018

    ##------ STEPS FOR HADD FOR THE TTBAR SAMPLES (OTHER SAMPLES WILL BE SIMILAR) -----##
    hadd -f TTToSemiLeptonic.root /hdfs/cms/store/user/rdewanje/tau_hlt_prod_2016/TTToSemiLeptonic*/*/*/000*/*.root
    mv TTToSemiLeptonic.root /hdfs/local/ram/run2_ntuples/full/2016/
    hadd -f TTTo2L2Nu.root /hdfs/cms/store/user/rdewanje/tau_hlt_prod_2016/TTTo2L2Nu*/*/*/000*/*.root
    mv TTTo2L2Nu.root /hdfs/local/ram/run2_ntuples/full/2016/

    hadd -f TTToSemiLeptonic.root /hdfs/cms/store/user/rdewanje/tau_hlt_prod_2017/TTToSemiLeptonic*/*/*/000*/*.root
    mv TTToSemiLeptonic.root /hdfs/local/ram/run2_ntuples/full/2017/
    hadd -f TTTo2L2Nu.root /hdfs/cms/store/user/rdewanje/tau_hlt_prod_2017/TTTo2L2Nu*/*/*/000*/*.root
    mv TTTo2L2Nu.root /hdfs/local/ram/run2_ntuples/full/2017/

    hadd -f TTToSemiLeptonic.root /hdfs/cms/store/user/rdewanje/tau_hlt_prod_2018/TTToSemiLeptonic*/*/*/000*/*.root
    mv TTToSemiLeptonic.root /hdfs/local/ram/run2_ntuples/full/2018/
    hadd -f TTTo2L2Nu.root /hdfs/cms/store/user/rdewanje/tau_hlt_prod_2018/TTTo2L2Nu*/*/*/000*/*.root
    mv TTTo2L2Nu.root /hdfs/local/ram/run2_ntuples/full/2018/
    .
    .
    .

    ##-------ALSO COPY THE PILEUP INFO FOR THE RESPECTIVE ERAS TO THIS DIRECTORY ------##
    cp TauTagAndProbe/data/2016/Pileup_Data2016.root /hdfs/local/ram/run2_ntuples
    cp TauTagAndProbe/data/2017/Pileup_Data2017.root /hdfs/local/ram/run2_ntuples
    cp TauTagAndProbe/data/2018/Pileup_Data2018.root /hdfs/local/ram/run2_ntuples


(4) After getting all the merged Run-2 Ntuples inside their respective (era-wise) sub-directories inside "run2_ntuples" folder execute this for skimming, Bkg. estimation, creation/fit of turn-ons and Tau Trigger SF computation:

    python processNtuples_stage2.py

(5a) The New Trigger scale factors (computed separately for true and fake taus) would be available inside the files named: "NewTriggerSFs_$ERA.root" (where, $ERA = 2016/2017/2018).

(5b) The Old Trigger scale factors (computed in 2 different ways: adding all fake-taus to DY MC process, subtracting all fake-taus from data) would be available inside the files named: "turn_on_$ERA_add-to-dy-mc_fitted_LATEST.root" and "turn_on_$ERA_subtract-from-data_fitted_LATEST.root" (where, $ERA = 2016/2017/2018).

(6) The plots for both old and new Trigger scale factors will be located inside the directory "Tau_Trigger_sf_plots" in .pdf, .root and .png formats for all Tau ID WPs, channels and tau decay modes.
