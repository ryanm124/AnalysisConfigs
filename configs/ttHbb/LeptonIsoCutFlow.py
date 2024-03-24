import matplotlib.pyplot as plt
from coffea.util import load
from pocket_coffea.utils.plot_utils import Shape
from IPython.display import display, HTML
display(HTML("<style>.container { width:95% !important; }</style>"))
display(HTML("<style>.output_png { height: auto; }</style>"))
import pickle
import numpy as np
import os
import awkward as ak
from omegaconf import OmegaConf
import math
import mplhep as hep
from coffea.nanoevents.methods import vector
from coffea.nanoevents.methods import candidate

o = load("/afs/cern.ch/user/a/asparker/public/ttHmain/PocketCoffea/AnalysisConfigs/configs/ttHbb/production_March20/output_March20_ALL_2018.coffea")
ElectronIsoWP = [0.40,0.35,0.30,0.25,0.20,0.15,0.10,0.05]
MuonPFIsoWP = [0.40,0.25,0.20,0.15,0.10,0.05]
MuonTkIsoWP = [0.10,0.05]
MuonMiniIsoWP = [0.40,0.20,0.10,0.05]
BBWorkingPoints = [0.0] #[0.3140,0.1566,0.0399] #Note, the newest paper from BTV on bbtagging references WPs but doesn't seem to list them
ElectronIsoCollections = ["jetRelIso", "pfRelIso03_all","miniPFRelIso_all", "pfRelIso03_chg","miniPFRelIso_chg", "dr03EcalRecHitSumEt", "dr03HcalDepth1TowerSumEt", "dr03TkSumPt", "dr03TkSumPtHEEP"]
MuonIsoCollections = ["jetRelIso","pfRelIso03_all", "miniPFRelIso_all" , "pfRelIso03_chg","miniPFRelIso_chg", "pfRelIso04_all"]
f = open("isoOptimizationTmux.txt", "w")

for muonIso in MuonIsoCollections:
    if "mini" in muonIso:
        muonWPs = MuonMiniIsoWP
    else:
        muonWPs = MuonPFIsoWP
    for eleIso in ElectronIsoCollections:
        p = {}
        parquetDir = "/afs/cern.ch/work/r/rmccarth/public/mergedFiles"
        for sample in os.listdir(parquetDir):
            if(("TTbb" in sample) or ("DATA" in sample)):
                continue
            p[sample] = {}
            if("v7" in sample):
                bbtag = "FatJetGood_btagDDBvL"
            else:
                bbtag = "FatJetGood_btagDDBvLV2"
            for cat in ["baseline"]:
                p[sample][cat] = ak.from_parquet(parquetDir+"/"+sample+"/"+cat+"/nominal/merge.parquet",columns=["ElectronGood_eta","ElectronGood_phi","ElectronGood_pt","ElectronGood_mass","ElectronGood_charge","MuonGood_eta","MuonGood_phi","MuonGood_pt","MuonGood_mass","MuonGood_charge","FatJetGood_eta","FatJetGood_phi","FatJetGood_pt","FatJetGood_mass","JetGood_eta","JetGood_phi","JetGood_pt",bbtag,"ElectronGood_"+eleIso,"MuonGood_"+muonIso,"weight_nominal"])
        for muonWP in muonWPs:
            for eleWP in ElectronIsoWP:
                for bbWP in BBWorkingPoints:
                    #p = {}
                    #parquetDir = "/eos/cms/store/user/asparker/ttHboosted_March20/mergedFiles"
                    cutFlow = {"Signal" : {"baseline": 0, "ee": 0, "emu": 0, "mumu": 0}, "TT" : {"baseline": 0, "ee": 0, "emu": 0, "mumu": 0}, "OtherBkg" : {"baseline": 0, "ee": 0, "emu": 0, "mumu": 0}}
                    for sample in os.listdir(parquetDir):
                        if(("TTbb" in sample) or ("DATA" in sample)):
                            continue
                        #print(sample,"muonIso",muonIso,"WP",muonWP,"eleIso",eleIso,"WP",eleWP,"bbWP",bbWP)
                        #p[sample] = {}
                        #for cat in os.listdir(parquetDir+"/"+sample):
                        for cat in ["baseline"]:
                            #p[sample][cat] = ak.from_parquet(parquetDir+"/"+sample+"/"+cat+"/nominal/merge.parquet")
                            eleMask = p[sample][cat]["ElectronGood_"+eleIso]<eleWP
                            muonMask = p[sample][cat]["MuonGood_"+muonIso]<muonWP
                            if("v7" in sample):
                                fatMask = p[sample][cat]["FatJetGood_btagDDBvL"]>bbWP
                            else:
                                fatMask = p[sample][cat]["FatJetGood_btagDDBvLV2"]>bbWP
                            lepEta = ak.concatenate((p[sample][cat]["ElectronGood_eta"][eleMask],p[sample][cat]["MuonGood_eta"][muonMask]),axis=1)
                            lepPhi = ak.concatenate((p[sample][cat]["ElectronGood_phi"][eleMask],p[sample][cat]["MuonGood_phi"][muonMask]),axis=1)
                            lepPt = ak.concatenate((p[sample][cat]["ElectronGood_pt"][eleMask],p[sample][cat]["MuonGood_pt"][muonMask]),axis=1)
                            lepMass = ak.concatenate((p[sample][cat]["ElectronGood_mass"][eleMask],p[sample][cat]["MuonGood_mass"][muonMask]),axis=1)
                            lepCharge = ak.concatenate((p[sample][cat]["ElectronGood_charge"][eleMask],p[sample][cat]["MuonGood_charge"][muonMask]),axis=1)
                            fatEta = p[sample][cat]["FatJetGood_eta"][fatMask]
                            fatPhi = p[sample][cat]["FatJetGood_phi"][fatMask]
                            fatPt = p[sample][cat]["FatJetGood_pt"][fatMask]
                            fatMass = p[sample][cat]["FatJetGood_mass"][fatMask]
                            jetEta = p[sample][cat]["JetGood_eta"]
                            jetPhi = p[sample][cat]["JetGood_phi"]
                            jetPt = p[sample][cat]["JetGood_pt"]
                            if(sample[-1]=="_"):
                                name = sample+ "2018"
                            else:
                                name = sample
                            weights = p[sample][cat]['weight_nominal'] / o['sum_genweights'][name]
                            
                            nlep = ak.num(lepEta,axis=1)
                            lep_mask = (nlep==2) & (ak.sum(lepPt>=25.0,axis=1)>=2)
                            lepEta = lepEta[lep_mask]
                            lepPhi = lepPhi[lep_mask]
                            lepPt = lepPt[lep_mask]
                            lepMass = lepMass[lep_mask]
                            lepCharge = lepCharge[lep_mask]
                            fatEta = fatEta[lep_mask]
                            fatPhi = fatPhi[lep_mask]
                            fatPt = fatPt[lep_mask]
                            fatMass = fatMass[lep_mask]
                            jetEta = jetEta[lep_mask]
                            jetPhi = jetPhi[lep_mask]
                            jetPt = jetPt[lep_mask]
                            weights = weights[lep_mask]

                            lep4Vec = ak.zip(
                            {
                            "pt": lepPt,
                            "eta": lepEta,
                            "phi": lepPhi,
                            "mass": lepMass,
                            "charge": lepCharge
                            },
                            with_name="PtEtaPhiMCandidate",
                            behavior= candidate.behavior,
                            )
                            ll = lep4Vec[:,0]+lep4Vec[:,1]
                            ll_charge = getattr(ll,"charge")
                            ll_mass = getattr(ll,"mass")
                            
                            min_mass = (ll_mass>20) & (ll_mass<76)
                            mass_window = ll_mass>106
                            mass_cut = (min_mass) | (mass_window)
                            dilep_mask = mass_cut & (ll_charge==0)
                            
                            lepEta = lepEta[dilep_mask]
                            lepPhi = lepPhi[dilep_mask]
                            lepPt = lepPt[dilep_mask]
                            lepMass = lepMass[dilep_mask]
                            lepCharge = lepCharge[dilep_mask]
                            fatEta = fatEta[dilep_mask]
                            fatPhi = fatPhi[dilep_mask]
                            fatPt = fatPt[dilep_mask]
                            fatMass = fatMass[dilep_mask]
                            jetEta = jetEta[dilep_mask]
                            jetPhi = jetPhi[dilep_mask]
                            jetPt = jetPt[dilep_mask]
                            weights = weights[dilep_mask]
                            
                            lep4Vec = ak.zip(
                            {
                            "pt": lepPt,
                            "eta": lepEta,
                            "phi": lepPhi,
                            "mass": lepMass,
                            },
                            with_name="PtEtaPhiMLorentzVector",
                            behavior= vector.behavior,
                            )
                            fat4Vec = ak.zip(
                            {
                            "pt": fatPt,
                            "eta": fatEta,
                            "phi": fatPhi,
                            "mass": fatMass,
                            },
                            with_name="PtEtaPhiMLorentzVector",
                            behavior= vector.behavior,
                            )
                            fatJetLeptonDeltaR = fat4Vec.metric_table(lep4Vec)
                            #print("lepEta",lepEta)
                            #print("fatEta",fatEta)
                            #print("fatjetLeptonEta",fatjetLeptonEta)
                            #print("fatjetLeptonPhi",fatjetLeptonPhi)
                            #print(fatJetLeptonDeltaR)
                            mask_lepton_cleaning = ak.prod(fatJetLeptonDeltaR > 0.8, axis=2) == 1
                            fatEta = fatEta[mask_lepton_cleaning]
                            fatPhi = fatPhi[mask_lepton_cleaning]
                            fatPt = fatPt[mask_lepton_cleaning]
                            fatMass = fatMass[mask_lepton_cleaning]
                            fatPt_mask = ak.sum(fatPt>=200,axis=1)>=1
                            
                            lepEta = lepEta[fatPt_mask]
                            lepPhi = lepPhi[fatPt_mask]
                            lepPt = lepPt[fatPt_mask]
                            lepMass = lepMass[fatPt_mask]
                            lepCharge = lepCharge[fatPt_mask]
                            fatEta = fatEta[fatPt_mask]
                            fatPhi = fatPhi[fatPt_mask]
                            fatPt = fatPt[fatPt_mask]
                            fatMass = fatMass[fatPt_mask]
                            jetEta = jetEta[fatPt_mask]
                            jetPhi = jetPhi[fatPt_mask]
                            jetPt = jetPt[fatPt_mask]
                            weights = weights[fatPt_mask]
                            #print(mask_lepton_cleaning)
                            lep4Vec = ak.zip(
                            {
                            "pt": lepPt,
                            "eta": lepEta,
                            "phi": lepPhi,
                            "mass": lepMass,
                            },
                            with_name="PtEtaPhiMLorentzVector",
                            behavior= vector.behavior,
                            )
                            jet4Vec = ak.zip(
                            {
                            "pt": jetPt,
                            "eta": jetEta,
                            "phi": jetPhi,
                            "mass": jetEta
                            },
                            with_name="PtEtaPhiMLorentzVector",
                            behavior= vector.behavior,
                            )
                            jetLeptonDeltaR = jet4Vec.metric_table(lep4Vec)
                            mask_lepton_cleaning = ak.prod(jetLeptonDeltaR > 0.4, axis=2) == 1
                            jetEta = jetEta[mask_lepton_cleaning]
                            jetPhi = jetPhi[mask_lepton_cleaning]
                            jetPt = jetPt[mask_lepton_cleaning]
                            fat4Vec = ak.zip(
                            {
                            "pt": fatPt,
                            "eta": fatEta,
                            "phi": fatPhi,
                            "mass": fatMass,
                            },
                            with_name="PtEtaPhiMLorentzVector",
                            behavior= vector.behavior,
                            )
                            jet4Vec = ak.zip(
                            {
                            "pt": jetPt,
                            "eta": jetEta,
                            "phi": jetPhi,
                            "mass": jetEta
                            },
                            with_name="PtEtaPhiMLorentzVector",
                            behavior= vector.behavior,
                            )
                            jetFatJetDeltaR = jet4Vec.metric_table(fat4Vec)
                            mask_fatjet_cleaning = ak.prod(jetFatJetDeltaR > 0.8, axis=2) == 1
                            #print("num jets removed",ak.sum(mask_fatjet_cleaning==0))
                            jetEta = jetEta[mask_fatjet_cleaning]
                            jetPhi = jetPhi[mask_fatjet_cleaning]
                            njet = ak.num(jetEta,axis=1)
                            
                            eventMask = (njet>=2)
                            #print("ll_charge",ll_charge)
                            #print("ll_mass",ll_mass)
                            #print("lepPt",lepPt)
                            #print("njet",njet)
                            #print("mass_cut",mass_cut)
                            #print("eventMask",eventMask)
                            #print("len eventMask",len(eventMask))
                            #print("num pass",ak.sum(eventMask))
                            
                            weights = weights[eventMask]
                            #print("weight sum",ak.sum(weights))
                            if("ttHTobb" in sample):
                                cutFlow["Signal"][cat] += ak.sum(weights)
                            elif("TTTo" in sample):
                                cutFlow["TT"][cat] += ak.sum(weights)
                            else:
                                cutFlow["OtherBkg"][cat] += ak.sum(weights)
                            
                    for cat in ["baseline"]:
                        f.write("muonIso: {} WP: {} eleIso: {} WP: {} bbWP: {} \n".format(muonIso,muonWP,eleIso,eleWP,bbWP))
                        f.write("Category: {}\n".format(cat))
                        signal = cutFlow["Signal"][cat]
                        bkg = cutFlow["OtherBkg"][cat]+cutFlow["TT"][cat]
                        if(bkg>0):
                            f.write("signal: {} TT: {} Other Bkg: {} TotalBkg: {} s/sqrt(b): {} s/b: {} \n".format(cutFlow["Signal"][cat],cutFlow["TT"][cat],cutFlow["OtherBkg"][cat],bkg,signal/math.sqrt(bkg),signal/bkg))
                            f.flush()
                        else:
                            f.write("signal: {} TT: {} Other Bkg: {} TotalBkg: {} \n".format(cutFlow["Signal"][cat],cutFlow["TT"][cat],cutFlow["OtherBkg"][cat],bkg))
                            f.flush()
f.close()
