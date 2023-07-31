from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel,get_nBtagEq,get_nBtagMin
from pocket_coffea.parameters.cuts import passthrough
from pocket_coffea.parameters.histograms import *
import os

from pocket_coffea.workflows.tthbb_base_processor import ttHbbBaseProcessor 

# importing custom cut functions
from custom_cut_functions import *
localdir = os.path.dirname(os.path.abspath(__file__))

# Loading default parameters
from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/params")

# merging additional analysis specific parameters
parameters = defaults.merge_parameters_from_files(default_parameters,
                                                  f"{localdir}/params/object_preselection.yaml",
                                                  f"{localdir}/params/triggers.yaml",
                                                  update=True)

# Configurator instance
cfg = Configurator(
    parameters = parameters,
    datasets = {
        "jsons": [f"{localdir}/datasets/backgrounds_MC_TTbb_dileptonic_redirector.json"
                    ],
        "filter" : {
            "samples": ["TTbbDiLeptonic"],
            "samples_exclude" : [],
            "year": ['2018']
        },
        "subsamples":{
            "TTToSemiLeptonic": {
                "=1b":  [get_nBtagEq(1, coll="Jet")],
                "=2b" : [get_nBtagEq(2, coll="Jet")],
                ">2b" : [get_nBtagMin(3, coll="Jet")]
            }
        }
    },

    workflow = ttHbbBaseProcessor,
    workflow_options = {},
    
    # Skimming and categorization
    skim = [
             get_nObj_min(4, 15., "Jet"),
             get_HLTsel()
             ],
             
    preselections = [dilepton_presel,
                      get_nObj_min(2,25,"LeptonGood")],
    
    categories = {
        "baseline": [passthrough],
        "1b" : [ get_nBtagEq(1, coll="BJetGood")],
        "2b" : [ get_nBtagEq(2, coll="BJetGood")],
        "3b" : [ get_nBtagEq(3, coll="BJetGood")],
        "4b" : [ get_nBtagEq(4, coll="BJetGood")]
    },
    
    # Weights configuration
    weights = {
        "common": {
            "inclusive": ["genWeight","lumi","XS",
                          "pileup",
                          "sf_ele_reco", "sf_ele_id",
                          "sf_mu_id","sf_mu_iso",
                          "sf_btag", "sf_jet_puId", 
                          ],
            "bycategory" : {
                "2jets_20pt" : []
            }
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


            }

        }

    },

    
    variables = {
        "HT" : HistConf([Axis(coll="events", field="events", bins=100, start=0, stop=200, label="HT")] ),
        "leading_jet_pt_eta" : HistConf(
            [
                Axis(coll="JetGood", field="pt", bins=40, start=0, stop=200, pos=0, label="Leading jet $p_T$"),
                Axis(coll="JetGood", field="eta", bins=40, start=-5, stop=5, pos=0, label="Leading jet $\eta$")
            ] ),
            
         #Plotting all jets together
         "all_jets_pt_eta" : HistConf(
            [
                Axis(coll="JetGood", field="pt", bins=40, start=0, stop=200, pos=None, label="All jets $p_T$"),
                Axis(coll="JetGood", field="eta", bins=40, start=-5, stop=5, pos=None, label="All jets $\eta$")
            ] ),
            
        "subleading_jetpt_MET" : HistConf(
            [
                Axis(coll="JetGood", field="pt", bins=40, start=0, stop=200, pos=0, label="Leading jet $p_T$"),
                Axis(coll="MET", field="pt", bins=40, start=0, stop=100, label="MET")
            ] ),
                
        **ele_hists(coll="ElectronGood", pos=0),
        **muon_hists(coll="MuonGood", pos=0),
        **count_hist(name="nElectronGood", coll="ElectronGood",bins=3, start=0, stop=3),
        **count_hist(name="nMuonGood", coll="MuonGood",bins=3, start=0, stop=3),
        **count_hist(name="nJets", coll="JetGood",bins=10, start=4, stop=14),
        **count_hist(name="nBJets", coll="BJetGood",bins=12, start=2, stop=14),
        **jet_hists(coll="JetGood", pos=0),
        **jet_hists(coll="JetGood", pos=1),
        **jet_hists(coll="JetGood", pos=2),
        
    },

    columns = {
        "common": {
             "inclusive": [],
             "bycategory": {}
        }

    }
)

run_options = {
        "executor"       : "dask/lxplus",
        "env"            : "singularity",
        "workers"        : 1,
        "scaleout"       : 50,
        "worker_image"   : "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-cc7-latest",
        "queue"          : "microcentury",
        "walltime"       : "00:40:00",
        "mem_per_worker" : "4GB", # GB
        "disk_per_worker" : "1GB", # GB
        "exclusive"      : False,
        "chunk"          : 400000,
        "retries"        : 50,
        "treereduction"  : 20,
        "adapt"          : False,
        
    }
