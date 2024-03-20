from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel,get_nBtagEq,get_nBtagMin, get_diLeptonFlavor
from pocket_coffea.parameters.cuts import passthrough
from pocket_coffea.parameters.histograms import *
from pocket_coffea.lib.columns_manager import ColOut
import os, ast
import workflow
from workflow import ttHbbBaseProcessor 

# importing custom cut functions
import custom_cut_functions
from  custom_cut_functions import *
import cloudpickle
cloudpickle.register_pickle_by_value(workflow)
cloudpickle.register_pickle_by_value(custom_cut_functions)
localdir = os.path.dirname(os.path.abspath(__file__))

from argParser import get_year_from_args
#import argparse
import os
from pocket_coffea.parameters import defaults

#parser = argparse.ArgumentParser(description='Description of your script')
#parser.add_argument('--year', type=int, help='Specify the year', required=True)
#args = parser.parse_args()

#aargs = get_year_from_args()
year = "2018" #c#aargs[0] #'2018'#args.year
samples = "bosons" #aargs[1]
samples_list_string = "WW,WZ,ZZ,WJetsToLNu_HT-200To400,WJetsToLNu_HT-400To600,WJetsToLNu_HT-600To800,WJetsToLNu_HT-800To1200,WJetsToLNu_HT-1200To2500,WJetsToLNu_HT-2500ToInf,WJetsToQQ_HT400to600_v7,WJetsToQQ_HT600to800_v7,WJetsToQQ_HT800toInf_v7,ZJetsToQQ_HT400to600_v7,ZJetsToQQ_HT600to800_v7,ZJetsToQQ_HT800toInf_v7,DYJetsToLL_M-50_v7" #aargs[2]

print("samples_list_string")
print(samples_list_string)
#samples_list = ast.literal_eval(samples_list_string)
samples_list = [sample.strip() for sample in samples_list_string.split(',')]
print(samples_list)


formatted_samples_files = []
formatted_samples = []
#formatted_samples_names = []

for i, s in enumerate(samples_list):
    formatted_samples_files.append(f"{localdir}/datasets/{s}_{year}.json")
    if 'v7' in s:
        formatted_samples.append(f"{s[:-3]}")
        print(s[:-3])
    else:
        formatted_samples.append(f"{s}")
        #formatted_samples_names.append(f"samples_{i+1}")

    
print("formatted_samples_files")
print(formatted_samples_files)

print("formatted_samples")
print(formatted_samples)

# Loading default parameters
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/params")

# merging additional analysis specific parameters
parameters = defaults.merge_parameters_from_files(default_parameters,
                                                  f"{localdir}/params/object_preselection.yaml",
                                                  f"{localdir}/params/triggers.yaml",
                                                  update=True)

# Configurator instance
datasets_jsons = [f"{formatted_sample}" for formatted_sample in formatted_samples_files]
datasets_filter_samples = [f"{formatted_sample_name}" for formatted_sample_name in formatted_samples ]



'''
formatted_samples_1 = None
formatted_samples_2 = None
formatted_samples_3 = None
formatted_samples_4 = None
formatted_samples_5 = None
print(samples)
for i,s in samples:
    if('_v7' in samples): s =s[:-3] # 'ZJetsToQQ_HT800toInf'
    print(s)
    if i == 0 :
        formatted_samples_1 = "{localdir}/datasets/{s}_{year}.json"
        samples_1 = s
    if i == 1 :
        formatted_samples_2 = "{localdir}/datasets/{s}_{year}.json"
      	samples_2 = s
        
formatted_samples = "{localdir}/datasets/{sample}_{year}.json"

#samples = aargs[1]
#year = '2018'
# Loading default parameters
from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/params")

# merging additional analysis specific parameters
parameters = defaults.merge_parameters_from_files(default_parameters,
                                                  f"{localdir}/params/object_preselection.yaml",
                                                  f"{localdir}/params/triggers.yaml",
                                                  update=True)
#sampleFilename = "DATA_SingleMuon_redirector.json"
#sampleName = "SingleMuon"
'''



# Configurator instance
#                  #f"{localdir}/datasets/DYJetsToLL_M-50_2018.json",#f"{localdir}/datasets/DATA_SingleMuon.json"                                                                               
cfg = Configurator(
    parameters = parameters,
    datasets = {
        "jsons": datasets_jsons ,
        "filter" : {
            "samples": datasets_filter_samples,
            "samples_exclude" : [],
            "year": [year]
        } 
    },

    workflow = ttHbbBaseProcessor,
    #workflow_options = {"dump_columns_as_arrays_per_chunk": "root://eosuser.cern.ch//eos/user/a/asparker/ttHbb/chunks"},
    workflow_options = {"dump_columns_as_arrays_per_chunk": "root://eoscms.cern.ch//eos/cms/store/user/asparker/ttHboosted/chunks"},
    
    # Skimming and categorization
    skim = [
             get_nObj_min(1, 200., "FatJet"),
             get_HLTsel(primaryDatasets=["DoubleEle","EleMu","DoubleMu"])
             ],
             
    preselections = [dilepton_presel,
                     get_nObj_min(2,25,"LeptonGood"),
                     get_nObj_min(1, 200., "FatJetGood")
                     ],
    
    categories = {
        "baseline": [passthrough],
        "ee" : [get_diLeptonFlavor("ee")],
        "emu" : [get_diLeptonFlavor("emu")],
        "mumu" : [get_diLeptonFlavor("mumu")]
    },
    
    # Weights configuration # "sf_ele_trigger", "sf_mu_trigger"
    weights = {
        "common": {
            "inclusive": ["genWeight","lumi","XS",
                          "pileup",
                          "sf_ele_reco", "sf_ele_id",
                          "sf_mu_id","sf_mu_iso",
                          "sf_btag", "sf_jet_puId" 
                          ],
         
        }
    },

    variations = {
        "weights": {
            "common": {
                "inclusive": [  "pileup",
                                "sf_ele_reco", "sf_ele_id",
                                "sf_mu_id", "sf_mu_iso",
                                 "sf_jet_puId","sf_btag"
                                
                              ]
            },
          
        },
       
    },

    
    variables = {
        **ele_hists(coll="ElectronGood", pos=0),
        **ele_hists(coll="ElectronGood", pos=1),
        **muon_hists(coll="MuonGood", pos=0),
        **muon_hists(coll="MuonGood", pos=1),
        **count_hist(name="nElectronGood", coll="ElectronGood",bins=3, start=0, stop=3),
        **count_hist(name="nMuonGood", coll="MuonGood",bins=3, start=0, stop=3),
        **count_hist(name="nJets", coll="JetGood",bins=8, start=0, stop=8),
        **count_hist(name="nBJets", coll="BJetGood",bins=8, start=0, stop=8),
        **count_hist(name="nBBFatJetGoodT", coll="BBFatJetGoodT",bins=3, start=0, stop=3),
        **count_hist(name="nBBFatJetGoodM", coll="BBFatJetGoodM",bins=3, start=0, stop=3),
        **count_hist(name="nBBFatJetGoodL", coll="BBFatJetGoodL",bins=3, start=0, stop=3),
        **jet_hists(coll="JetGood", pos=0),
        **jet_hists(coll="JetGood", pos=1),
        **jet_hists(coll="JetGood", pos=2),
        **jet_hists(coll="JetGood", pos=3),
        **jet_hists(coll="JetGood", pos=4),
        **jet_hists(name="bjet",coll="BJetGood", pos=0),
        **jet_hists(name="bjet",coll="BJetGood", pos=1),
        **jet_hists(name="bjet",coll="BJetGood", pos=2),
        **fatjet_hists(name="fatjet",coll="FatJetGood"),
        **fatjet_hists(name="bbfatjetTight",coll="BBFatJetGoodT"),
        **fatjet_hists(name="bbfatjetMedium",coll="BBFatJetGoodM"),
        **fatjet_hists(name="bbfatjetLoose",coll="BBFatJetGoodL"),

       # 2D plots
       "jet_eta_pt_leading": HistConf(
           [
               Axis(coll="JetGood", field="pt", pos=0, bins=40, start=0, stop=1000,
                    label="Leading jet $p_T$"),
               Axis(coll="JetGood", field="eta", pos=0, bins=40, start=-2.4, stop=2.4,
                    label="Leading jet $\eta$"),
           ]
       ),
       "jet_eta_pt_all": HistConf(
           [
               Axis(coll="JetGood", field="pt", bins=40, start=0, stop=1000,
                    label="jet $p_T$"),
               Axis(coll="JetGood", field="eta", bins=40, start=-2.4, stop=2.4,
                    label="jet $\eta$")
           ]
       ),

        "subleading_jetpt_MET" : HistConf(
            [
                Axis(coll="JetGood", field="pt", bins=40, start=0, stop=200, pos=1, label="Subleading jet $p_T$"),
                Axis(coll="MET", field="pt", bins=40, start=0, stop=100, label="MET")
            ] ),

    },
    
    columns = {
        "common": {
            "inclusive": [
                ColOut("JetGood", ["eta","pt","phi","btagDeepFlavB"]),
                ColOut("FatJetGood", ["eta", "pt", "phi", "mass", "msoftdrop", "tau1", "tau2", "tau3", "tau4", "btagDDBvLV2", "deepTagMD_ZHbbvsQCD", "deepTagMD_ZHccvsQCD", "deepTagMD_HbbvsQCD", "deepTagMD_bbvsLight", "btagHbb", "rhoQCD"]),
                ColOut("MuonGood",["eta","pt","phi","mass","charge","jetRelIso","pfRelIso03_all", "miniPFRelIso_all" , "pfRelIso03_chg","miniPFRelIso_chg", "mvaTTH", "pfRelIso04_all" ]),
                ColOut("ElectronGood",["eta","pt","phi","mass","charge","jetRelIso", "pfRelIso03_all","miniPFRelIso_all", "pfRelIso03_chg","miniPFRelIso_chg", "mvaTTH","dr03EcalRecHitSumEt", "dr03HcalDepth1TowerSumEt", "dr03TkSumPt", "dr03TkSumPtHEEP" ]),                
                ColOut("events",["genTtbarId"],store_size=False),
                ColOut("MET",["pt"])
            ]
        }
    }
    
    
 
 
 
 
 
 
 
 
 
)
'''
run_options = {
        "executor"       : "dask/lxplus",
        "env"            : "singularity",
        "workers"        : 1,
        "scaleout"       : 50,
        "worker_image"   : "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-cc7-latest",
        "queue"          : "workday", # workday = 8 hour max job runtime : https://batchdocs.web.cern.ch/local/submit.html
        "walltime"       : "00:40:00",
        "mem_per_worker" : "4GB", # GB
        "disk_per_worker" : "1GB", # GB
        "exclusive"      : False,
        "chunk"          : 400000,
        "retries"        : 50,
        "treereduction"  : 20,
        "adapt"          : False,
        "skipbadfiles"   : 10        
    }
'''
