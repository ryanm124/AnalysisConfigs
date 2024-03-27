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
import mplhep as hep
from coffea.nanoevents.methods import vector
from coffea.nanoevents.methods import candidate
#/afs/cern.ch/user/a/asparker/public/ttHmain
parameters_dump = OmegaConf.load("/afs/cern.ch/user/a/asparker/public/ttHmain/PocketCoffea/AnalysisConfigs/configs/ttHbb/params/plotting_style.yaml")

def getName(sam_name):
    sam_namey = sam_name+'__2018'
    if 'ttHToNonbb' in sam_name: sam_namey ='ttHToNonbb_2018'
    if 'DYJetsToLL_M-50' in sam_name: sam_namey ='DYJetsToLL_M-50_v7__2018'
    if 'ZJetsToQQ' in sam_name: sam_namey = sam_name + '_v7__2018'
    if 'ST' in sam_name: sam_namey = sam_name + '_v7__2018'
    if 'QCD' in sam_name: sam_namey = sam_name + '_v7__2018'
    if 'WJets' in sam_name: sam_namey = sam_name + '_v7__2018'
    if 'WJetsToLNu' in sam_name: sam_namey = sam_name + '__2018'
    if 'TTWJets' in sam_name: sam_namey = sam_name + '__2018'
    if 'TTGJets' in sam_name: sam_namey = sam_name + '__2018'
    if 'THW' in sam_name: sam_namey = sam_name + '__2018'
    if 'WW' in sam_name: sam_namey = sam_name + '__2018'
    if 'WZ' in sam_name: sam_namey = sam_name + '__2018'
    if 'ZZ' in sam_name: sam_namey = sam_name + '__2018'
    if 'ttHTobb' in sam_name: sam_namey = sam_name + '_2018'
    return sam_namey

def getPlotParams(var):
    #return format: bins, xmin, xmax, label
    col = var.split("_")[0]
    quantity = var.split("_")[1]
    if quantity == 'eta':
        return 30, -2.5, 2.5, fr"{col} $\eta$"
    if quantity == 'phi':
        return 32, -math.pi, math.pi, fr"{col} $\phi$"
    if ("FatJet" in col) and not ("N" in quantity) and not ("Flavour" in quantity):
        if quantity == 'pt':
            return 100, 0, 1000, r"FatJet $p_{T}$ [GeV]"
        elif quantity == 'mass':
            return 50, 0, 400, r"FatJet mass [GeV]"
        elif quantity == 'msoftdrop':
            return 50, 0, 400, r"FatJet $m_{SD}$ [GeV]"
        elif quantity == 'rho':
            return 100, -8, 0, r"FatJet $\rho$"
        elif quantity == 'rhoQCD':
            return 100, -3, 0.5, r"FatJet $\rho$QCD"
        elif quantity == 'tau1':
            return 20, 0, 1, r"$\tau_{1}$"
        elif quantity == 'tau2':
            return 20, 0, 1, r"$\tau_{2}$"
        elif quantity == 'tau3':
            return 20, 0, 1, r"$\tau_{3}$"
        elif quantity == 'tau4':
            return 20, 0, 1, r"$\tau_{4}$"
        elif quantity == 'btagDDBvLV2':
            return 50, 0, 1, f"{var}"
        else:
            return -1, -1, -1, f"{var}"
    elif ("Jet" in col) and not ("Flavour" in quantity):
        if quantity == 'pt':
            return 100, 0, 1000, r"Jet $p_{T}$ [GeV]"
        elif quantity == 'btagDeepFlavB':
            return 30, 0, 1, "AK4 DeepJet b-tag score"
        elif quantity == 'mass':
            return 50, 0, 400, r"Jet mass [GeV]"
        elif quantity == 'N':
            return 8, 0, 8, fr"{col} N"
        else:
            return -1, -1, -1, f"{var}"
    elif ("Lepton" in col):
        if quantity == 'pt':
            return 30, 0, 300, r"Lepton $p_{T}$"
        else:
            return -1, -1, -1, f"{var}"
    else:
        return -1, -1, -1, f"{var}"


def deltaR(phis,etas):
    leftPhi,rightPhi = ak.unzip(phis)
    leftEta,rightEta = ak.unzip(etas)
    return np.sqrt(pow(leftEta-rightEta,2)+pow(leftPhi-rightPhi,2))

o = load("/afs/cern.ch/user/a/asparker/public/ttHmain/PocketCoffea/AnalysisConfigs/configs/ttHbb/production_March20/output_March20_ALL_2018.coffea")
p = {}
parquetDir = "/eos/cms/store/user/asparker/ttHboosted_March20/mergedFiles"
for sample in os.listdir(parquetDir):
    p[sample] = {}
    #for cat in os.listdir(parquetDir+"/"+sample):
    for cat in ["baseline"]:
        p[sample][cat] = ak.from_parquet(parquetDir+"/"+sample+"/"+cat+"/nominal/merge.parquet")
        
ElectronIsoWP = [0.06] 
MuonPFIsoWP = [0.25] 
MuonTkIsoWP = [0.10,0.05]
MuonMiniIsoWP = [0.40,0.20,0.10,0.05]
BBWorkingPoints = [0.0] 
ElectronIsoCollections = ["pfRelIso03_all" ]
MuonIsoCollections =  ["pfRelIso04_all"] 
f = open("isoOptimizationFast_oldcuts.txt", "a")

for muonIso in MuonIsoCollections:
    if "mini" in muonIso:
        muonWPs = MuonMiniIsoWP
    else:
        muonWPs = MuonPFIsoWP
    for eleIso in ElectronIsoCollections:
        for muonWP in muonWPs:
            for eleWP in ElectronIsoWP:
                for bbWP in BBWorkingPoints:
                    #p = {}
                    parquetDir = "/eos/cms/store/user/asparker/ttHboosted_March20/mergedFiles"
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
                            jetEta = p[sample][cat]["JetGood_eta"]
                            jetPhi = p[sample][cat]["JetGood_phi"]
                            #print("len jetPhi",len(jetPhi))
                            fatjetLeptonEta = ak.cartesian([fatEta,lepEta],axis=1,nested=True)
                            fatjetLeptonPhi = ak.cartesian([fatPhi,lepPhi],axis=1,nested=True)
                            fatJetLeptonDeltaR = deltaR(fatjetLeptonPhi,fatjetLeptonEta)
                            #print("lepEta",lepEta)
                            #print("fatEta",fatEta)
                            #print("fatjetLeptonEta",fatjetLeptonEta)
                            #print("fatjetLeptonPhi",fatjetLeptonPhi)
                            #print(fatJetLeptonDeltaR)
                            mask_lepton_cleaning = ak.prod(fatJetLeptonDeltaR > 0.8, axis=2) == 1
                            fatEta = fatEta[mask_lepton_cleaning]
                            fatPhi = fatPhi[mask_lepton_cleaning]
                            fatPt = fatPt[mask_lepton_cleaning]
                            #print(mask_lepton_cleaning)
                            jetLeptonEta = ak.cartesian([jetEta,lepEta],axis=1,nested=True)
                            jetLeptonPhi = ak.cartesian([jetPhi,lepPhi],axis=1,nested=True)
                            jetLeptonDeltaR = deltaR(jetLeptonPhi,jetLeptonEta)
                            mask_lepton_cleaning = ak.prod(fatJetLeptonDeltaR > 0.4, axis=2) == 1
                            jetEta = jetEta[mask_lepton_cleaning]
                            jetPhi = jetPhi[mask_lepton_cleaning]
                            jetFatJetEta = ak.cartesian([jetEta,fatEta],axis=1,nested=True)
                            jetFatJetPhi = ak.cartesian([jetPhi,fatPhi],axis=1,nested=True)
                            jetFatJetDeltaR = deltaR(jetFatJetPhi,jetFatJetEta)
                            mask_fatjet_cleaning = ak.prod(jetFatJetDeltaR > 0.8, axis=2) == 1
                            #print("num jets removed",ak.sum(mask_fatjet_cleaning==0))
                            jetEta = jetEta[mask_fatjet_cleaning]
                            jetPhi = jetPhi[mask_fatjet_cleaning]
                            #print("len jetPhi end",len(jetPhi))
                            lepEta = ak.pad_none(lepEta,2)
                            lepPhi = ak.pad_none(lepPhi,2)
                            lepPt = ak.pad_none(lepPt,2)
                            lepMass = ak.pad_none(lepMass,2)
                            lepCharge = ak.pad_none(lepCharge,2)
                            nlep = ak.num(lepEta[~ak.is_none(lepEta,axis=1)])
                            #print("nlep",nlep)
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
                            ll_charge = ak.where((nlep==2),getattr(ll,"charge"),-1)
                            ll_mass = ak.where((nlep==2),getattr(ll,"mass"),-1)
                            njet = ak.num(jetEta,axis=1)
                            min_mass = (ll_mass>20) & (ll_mass<76)
                            mass_window = ll_mass>106
                            mass_cut = (min_mass) | (mass_window)
                            eventMask = (nlep==2) & (ak.sum(lepPt>=25.0,axis=1)>=2) & (njet>=2) & (ll_charge==0) & (mass_cut) & (ak.sum(fatPt>=200,axis=1)>=1)
                            #print("ll_charge",ll_charge)
                            #print("ll_mass",ll_mass)
                            #print("lepPt",lepPt)
                            #print("njet",njet)
                            #print("mass_cut",mass_cut)
                            #print("eventMask",eventMask)
                            #print("len eventMask",len(eventMask))
                            #print("num pass",ak.sum(eventMask))
                            if(sample[-1]=="_"):
                                name = sample+ "2018"
                            else:
                                name = sample
                            weights = p[sample][cat]['weight_nominal'] / o['sum_genweights'][name]
                            weights = weights[eventMask]
                            print("weight sum",ak.sum(weights))
                            if("ttHTobb" in sample):
                                cutFlow["Signal"][cat] += ak.sum(weights)
                            elif("TTTo" in sample):
                                cutFlow["TT"][cat] += ak.sum(weights)
                            else:
                                cutFlow["OtherBkg"][cat] += ak.sum(weights)
                    for cat in ["baseline"]:
                        print("muonIso",muonIso,"WP",muonWP,"eleIso",eleIso,"WP",eleWP,"bbWP",bbWP)
                        f.write("muonIso: {} WP: {} eleIso: {} WP: {} bbWP: {} \n".format(muonIso,muonWP,eleIso,eleWP,bbWP))
                        print("Category:",cat)
                        f.write("Category: {}\n".format(cat))
                        signal = cutFlow["Signal"][cat]
                        bkg = cutFlow["OtherBkg"][cat]+cutFlow["TT"][cat]
                        if(bkg>0):
                            print("signal:",cutFlow["Signal"][cat],"TT:",cutFlow["TT"][cat],"Other Bkg",cutFlow["OtherBkg"][cat],"TotalBkg",bkg,"s/sqrt(b)",signal/math.sqrt(bkg),"s/b",signal/bkg)
                            f.write("signal: {} TT: {} Other Bkg: {} TotalBkg: {} s/sqrt(b): {} s/b: {} \n".format(cutFlow["Signal"][cat],cutFlow["TT"][cat],cutFlow["OtherBkg"][cat],bkg,signal/math.sqrt(bkg),signal/bkg))
                            f.flush()
                        else:
                            print("signal:",cutFlow["Signal"][cat],"TT:",cutFlow["TT"][cat],"Other Bkg",cutFlow["OtherBkg"][cat],"TotalBkg",bkg)
                            f.write("signal: {} TT: {} Other Bkg: {} TotalBkg: {} \n".format(cutFlow["Signal"][cat],cutFlow["TT"][cat],cutFlow["OtherBkg"][cat],bkg))
                            f.flush()

f.close()
