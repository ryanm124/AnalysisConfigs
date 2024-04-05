


import matplotlib.pyplot as plt
from coffea.util import load
from pocket_coffea.utils.plot_utils import Shape
from IPython.display import display, HTML
import pickle
import numpy as np
import os
import awkward as ak
from omegaconf import OmegaConf
import math
import pyarrow as pa

def deltaR(phis,etas):
    leftPhi,rightPhi = ak.unzip(phis)
    leftEta,rightEta = ak.unzip(etas)
    return np.sqrt(pow(leftEta-rightEta,2)+pow(leftPhi-rightPhi,2))


samples=(
    "DATA_SingleMuon_2018_EraA",
    "DATA_SingleMuon_2018_EraB",
    "DATA_SingleMuon_2018_EraC",
    "DATA_SingleMuon_2018_EraD",
    "DATA_EGamma_2018_EraA",
    "DATA_EGamma_2018_EraB",
    "DATA_EGamma_2018_EraC",
    "DATA_EGamma_2018_EraD",
)
'''
    "DATA_SingleMuon_2018_EraA",
    "DATA_SingleMuon_2018_EraB",
    "DATA_SingleMuon_2018_EraC",
    "DATA_SingleMuon_2018_EraD"
)    
    "QCD_HT500to700_v7",
         "QCD_HT700to1000_v7",
         "QCD_HT1000to1500_v7",
         "QCD_HT1500to2000_v7",
         "QCD_HT2000toInf_v7"
        )
'''

years = ["2018"]
cats = ["baseline","ee","emu","mumu"]
variations = ["nominal"]
#filepath = "/eos/cms/store/user/asparker/ttHboosted_March20_2/chunks/"
#filepath="/eos/user/r/rmccarth/ttHbb/chunks/"
#outputpath = "/eos/cms/store/user/asparker/ttHboosted_March20/mergedFiles/"
filepath = "/eos/cms/store/user/asparker/ttHboosted_april3/chunks/"
outputpath = "/afs/cern.ch/work/a/asparker/public/ttHboosted_april3/chunks/"


year = "2018"
for sample in samples: #["ZZ__2018","WJetsToLNu_HT-200To400__2018","WJetsToLNu_HT-2500ToInf__2018","ST_tW_antitop_5f_inclusiveDecays_v7__2018","QCD_HT1500to2000_v7__2018","QCD_HT2000toInf_v7__2018","TTToSemiLeptonic__2018"]:
    #if ("ttHTo" or "Era" ) in sample :sample+="_"
    #else:
    #sample+="_"
    for cat in cats:
        for variation in variations:
            path = filepath+sample+"/"+cat+"/"+variation
            print(path)
            outpath = outputpath+sample+"/"+cat+"/"+variation
            if(os.path.isdir(outpath)):
                continue
            ak.to_parquet.dataset(path)
            output = ak.from_parquet(path)

            output['FatJetGood_tau21'] = np.divide(output['FatJetGood_tau2'],output['FatJetGood_tau1'])
            output['FatJetGood_tau32'] = np.divide(output['FatJetGood_tau3'],output['FatJetGood_tau2'])
            output['FatJetGood_tau31'] = np.divide(output['FatJetGood_tau3'],output['FatJetGood_tau1'])
            output['FatJetGood_tau43'] = np.divide(output['FatJetGood_tau4'],output['FatJetGood_tau3'])

            '''
            output['BBFatJetGoodL_tau21'] = np.divide(output['BBFatJetGoodL_tau2'],output['BBFatJetGoodL_tau1'])
            output['BBFatJetGoodL_tau32'] = np.divide(output['BBFatJetGoodL_tau3'],output['BBFatJetGoodL_tau2'])
            output['BBFatJetGoodL_tau31'] = np.divide(output['BBFatJetGoodL_tau3'],output['BBFatJetGoodL_tau1'])
            output['BBFatJetGoodL_tau43'] = np.divide(output['BBFatJetGoodL_tau4'],output['BBFatJetGoodL_tau3'])
                
            output['BBFatJetGoodM_tau21'] = np.divide(output['BBFatJetGoodM_tau2'],output['BBFatJetGoodM_tau1'])
            output['BBFatJetGoodM_tau32'] = np.divide(output['BBFatJetGoodM_tau3'],output['BBFatJetGoodM_tau2'])
            output['BBFatJetGoodM_tau31'] = np.divide(output['BBFatJetGoodM_tau3'],output['BBFatJetGoodM_tau1'])
            output['BBFatJetGoodM_tau43'] = np.divide(output['BBFatJetGoodM_tau4'],output['BBFatJetGoodM_tau3'])

            output['BBFatJetGoodT_tau21'] = np.divide(output['BBFatJetGoodT_tau2'],output['BBFatJetGoodT_tau1'])
            output['BBFatJetGoodT_tau32'] = np.divide(output['BBFatJetGoodT_tau3'],output['BBFatJetGoodT_tau2'])
            output['BBFatJetGoodT_tau31'] = np.divide(output['BBFatJetGoodT_tau3'],output['BBFatJetGoodT_tau1'])
            output['BBFatJetGoodT_tau43'] = np.divide(output['BBFatJetGoodT_tau4'],output['BBFatJetGoodT_tau3'])

            jetLeptonEta = ak.cartesian([output["LeptonGood_eta"],output["JetGood_eta"]],axis=1)
            jetLeptonPhi = ak.cartesian([output["LeptonGood_phi"],output["JetGood_phi"]],axis=1)
            jetFatJetEta = ak.cartesian([output["FatJetGood_eta"],output["JetGood_eta"]],axis=1)
            jetFatJetPhi = ak.cartesian([output["FatJetGood_phi"],output["JetGood_phi"]],axis=1)
            fatjetLeptonEta = ak.cartesian([output["LeptonGood_eta"],output["FatJetGood_eta"]],axis=1)
            fatjetLeptonPhi = ak.cartesian([output["LeptonGood_phi"],output["FatJetGood_phi"]],axis=1)
            leptonsPhi = ak.to_list(ak.combinations(output["LeptonGood_phi"],2))
            leptonsEta = ak.to_list(ak.combinations(output["LeptonGood_eta"],2))
            jetPhi = ak.to_list(ak.combinations(output["JetGood_phi"],2))
            jetEta = ak.to_list(ak.combinations(output["JetGood_eta"],2))

            if(len(leptonsPhi)!=0 and len(leptonsEta)!=0):
                leptonsDeltaR = deltaR(leptonsPhi,leptonsEta)
            else:
                leptonsDeltaR = []
            output['LeptonGoodCombinations_deltaR'] = leptonsDeltaR

            if(len(jetPhi)!=0 and len(jetEta)!=0):
                jetDeltaR = deltaR(jetPhi,jetEta)
            else:
                jetDeltaR = []
            output['JetGoodCombinations_deltaR'] = jetDeltaR

            if(len(leptonsPhi)!=0 and len(leptonsEta)!=0 and len(jetPhi)!=0 and len(jetEta)!=0):
                jetLeptonDeltaR = deltaR(jetLeptonPhi,jetLeptonEta)
            else:
                jetLeptonDeltaR = []
            output['JetLeptonGoodCombinations_deltaR'] = jetLeptonDeltaR

            if(len(fatjetLeptonPhi)!=0 and len(fatjetLeptonEta)!=0):
                fatjetLeptonDeltaR = deltaR(fatjetLeptonPhi,fatjetLeptonEta)
            else:
                fatjetLeptonDeltaR = []
            output['FatJetLeptonGoodCombinations_deltaR'] = fatjetLeptonDeltaR

            if(len(jetFatJetPhi)!=0 and len(jetFatJetEta)!= 0):
                jetFatJetDeltaR = deltaR(jetFatJetPhi,jetFatJetEta)
            else:
                jetFatJetDeltaR = []
            output['JetFatJetGoodCombinations_deltaR'] = jetFatJetDeltaR
            '''
            os.makedirs(outpath)
            ak.to_parquet(output,f"{outpath}/merge.parquet")


            
