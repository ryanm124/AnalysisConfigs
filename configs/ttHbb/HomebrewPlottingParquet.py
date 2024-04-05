
import matplotlib.pyplot as plt
from coffea.util import load
from pocket_coffea.utils.plot_utils import Shape
from IPython.display import display, HTML
display(HTML("<style>.container { width:95% !important; }</style>"))
display(HTML("<style>.output_png { height: auto; }</style>"))
import pickle
import numpy as np
import os
if not os.path.exists('hists'):
    os.makedirs('hists')
import awkward as ak
from omegaconf import OmegaConf
import math
import mplhep as hep
from coffea.nanoevents.methods import vector
from coffea.nanoevents.methods import candidate

parameters_dump = OmegaConf.load("/afs/cern.ch/work/a/asparker/public/ttHwork/PocketCoffea/AnalysisConfigs/configs/ttHbb/params/plotting_style.yaml")

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



o = load("/afs/cern.ch/work/a/asparker/public/ttHwork/PocketCoffea/AnalysisConfigs/configs/ttHbb/april3/output_most_2018.coffea")
#("/afs/cern.ch/user/a/asparker/public/ttHmain/PocketCoffea/AnalysisConfigs/configs/ttHbb/production_March20/output_March20_ALL_2018.coffea")  

print(o['sum_genweights'].keys())

p = {}
parquetDir ="/afs/cern.ch/work/a/asparker/public/ttHboosted_april3/chunks"

#o = load("/afs/cern.ch/user/a/asparker/public/ttHmain/PocketCoffea/AnalysisConfigs/configs/ttHbb/production_March20/output_March20_ALL_2018.coffea")

#print(o['sum_genweights'].keys())

#p = {}
#parquetDir = "/eos/cms/store/user/asparker/ttHboosted_april3/chunks"
#"/afs/cern.ch/work/a/asparker/public/ttHboosted_March20_merge" # "root://eoscms.cern.ch//eos/cms/store/user/asparker/ttHboosted_March20/mergedFiles" # 
for sample in os.listdir(parquetDir):
#for sample in ['TTbbDiLeptonic_Powheg__']:

    print(sample)
    p[sample] = {}
    #for cat in os.listdir(parquetDir+"/"+sample):
    for cat in ["baseline"]:
        p[sample][cat] = ak.from_parquet(parquetDir+"/"+sample+"/"+cat+"/nominal/merge.parquet")

print(p.keys())

#User parameters
#cats = ["baseline","ee","emu","mumu"]
cats = ["baseline"]
colors = ["red","orange","gold","yellow","green","cyan","blue","purple","sienna","greenyellow","lightblue"]
vars = ["JetGood_pt"]#[] #columns to plot, leave empty for all columns
modifier = "pt_1" #which objects to plot per event, options are: 'pt_#' where # goes from 1 to N with 1 being the leading pt object
                  #                                              'max' which plots the highest value of the parameter per event
                  #                                              'all' which plots all values of the parameter per event

dataSamples = []
for sample in p.keys():
    if 'DATA' in sample:
        dataSamples.append(sample)

for cat in cats:
    if(not len(vars)):
        vars = p['ttHTobb_']['baseline'].fields
    for var in vars:
        if(("weight" in var) or ("Idx" in var) or ("Ttbar" in var)):
            continue
        print("2018", cat, var)
        col = var.split("_")[0]
        quantity = var.split("_")[1]
        all_MC = {}
        all_weight = {}
        for label, samples in parameters_dump['plotting_style']['samples_groups'].items():
            label = parameters_dump['plotting_style']['labels_mc'].get(label, label)
            label_data = []
            label_weight = {}
            for sample in samples:
                print(sample)
                
                weight = {}
                ttbb_weight = {}
                name = getName(sample)
                if(("V2" in var) and ("v7" in name)):
                    data = p[name][cat][var[:-2]]
                else:
                    data = p[name[:-4]][cat][var]
                    #name = name[:-4]
                metpt = p[name[:-4]][cat]["MET_pt"] 
                met_mask = metpt > 20
                #data = data[~met_mask]    
                for weightVar in p[name[:-4]][cat].fields:          
                    if ("weight" in weightVar):
                        weight[weightVar] = p[name[:-4]][cat][weightVar]#[~met_mask]
                        weight[weightVar] = weight[weightVar] / o['sum_genweights'][name]
                if(sample=="TTTo2L2Nu" or sample=="TTToSemiLeptonic" or sample=="TTToHadronic"):
                    if sample=="TTTo2L2Nu":
                        ttbbSample = 'TTbbDiLeptonic_Powheg__'
                    elif sample=="TTToSemiLeptonic":
                        ttbbSample = 'TTbbSemiLeptonic_Powheg__'
                    else:
                        ttbbSample = 'TTbbHadronic_Powheg__'
                    genTtbarId = p[name[:-4]][cat]["events_genTtbarId"] % 100
                    tt_ttb_mask = genTtbarId > 50
                    
                    if(("V2" in var) and ("v7" in name)):
                        ttbb_data = p[ttbbSample][cat][var[:-2]]
                    else:
                        ttbb_data = p[ttbbSample][cat][var]
                    genTtbarId = p[ttbbSample][cat]["events_genTtbarId"] % 100
                    ttbb_ttb_mask = genTtbarId > 50
                    for weightVar in p[name[:-4]][cat].fields:
                        if ("weight" in weightVar):
                            ttbb_weight[weightVar] = p[ttbbSample][cat][weightVar]#[~met_mask]
                            print("....")
                            ttbb_weight[weightVar] = ttbb_weight[weightVar] / o['sum_genweights'][ttbbSample+'2018']
                            ttbb_weight[weightVar] = ttbb_weight[weightVar][ttbb_ttb_mask]
                            C = ak.sum(weight[weightVar][tt_ttb_mask])
                            D = ak.sum(ttbb_weight[weightVar])
                            ttbb_weight[weightVar] = ttbb_weight[weightVar] * (C/D)

                    for key in weight:
                        weight[key] = weight[key][~tt_ttb_mask]#[~met_mask]     
                    
                    if(quantity!="N" and col!="events"):
                        print(name[:-4])
                        pt_data = p[name[:-4]][cat][col+"_pt"]#[~met_mask]
                        sortIndices = ak.argsort(pt_data,ascending=False)
                        data = data[sortIndices]
                        ttbb_pt_data = p[ttbbSample][cat][col+"_pt"]#[~met_mask]
                        sortIndices = ak.argsort(ttbb_pt_data,ascending=False)
                        ttbb_data = ttbb_data[sortIndices]
                        if "pt" in modifier:
                            index = int(modifier.split("_")[1]) - 1
                            ttbb_data = ak.mask(ttbb_data, ak.num(ttbb_data) > index)[:, index]
                            data = ak.mask(data, ak.num(data) > index)[:, index]
                        elif modifier=="max":
                            data = ak.max(data,axis=1)
                            ttbb_data = ak.max(ttbb_data,axis=1)
                    data = data[~tt_ttb_mask]#[~met_mask]
                    ttbb_data = ttbb_data[ttbb_ttb_mask]#[~met_mask]
                    data = ak.concatenate((data,ttbb_data),axis=0)
                    for key in weight:
                        weight[key] = ak.concatenate((weight[key],ttbb_weight[key]),axis=0)
                        weight[key] = weight[key]#[~met_mask]
                elif(quantity!="N" and col!="events"):
                    pt_data = p[name[:-4]][cat][col+"_pt"]#[~met_mask]
                    sortIndices = ak.argsort(pt_data,ascending=False)
                    data = data[sortIndices]
                    if "pt" in modifier:
                        index = int(modifier.split("_")[1]) - 1
                        data = ak.mask(data, ak.num(data) > index)[:, index]
                    elif modifier=="max":
                        data = ak.max(data,axis=1)
                for key in weight:
                    weight[key], data = ak.broadcast_arrays(weight[key],data)
                    if(modifier=="all" and quantity!="N" and col!="events"):
                        weight[key] = ak.flatten(weight[key])
                if(modifier=="all" and quantity!="N" and col!="events"): 
                    data = ak.flatten(data)
                none_mask = ak.is_none(data)
                data = data[~none_mask]
                for key in weight:
                    weight[key] = weight[key][~none_mask]
                inf_mask = np.isinf(data)
                data = data[~inf_mask]
                for key in weight:
                    weight[key] = weight[key][~inf_mask]
                label_data = ak.concatenate((label_data,data),axis=0)
                for key in weight:
                    if key not in label_weight:
                        label_weight[key] = []
                    label_weight[key] = ak.concatenate((label_weight[key],weight[key]),axis=0)
            all_MC[label] = label_data
            for key in label_weight:
                if key not in all_weight:
                    all_weight[key] = {}
                all_weight[key][label] = label_weight[key]

        bins, xmin, xmax, xlabel = getPlotParams(var)
        plt.style.use([hep.style.ROOT, {'font.size': 22}])
        f2, a2 = plt.subplots()
        if not ("Gen" in var) and not ("gen" in var):
            fig, (ax, rax) = plt.subplots(2, 1, figsize=[12,12], gridspec_kw={"height_ratios": [3,1]}, sharex=True)
            fig.subplots_adjust(hspace=0.06)
            if(bins==-1):
                n, bin_edges, patches = ax.hist(list(all_MC.values()),weights=list(all_weight["weight_nominal"].values()),stacked=True,label=list(all_MC.keys()),color=colors)
            else:
                n, bin_edges, patches = ax.hist(list(all_MC.values()),weights=list(all_weight["weight_nominal"].values()),stacked=True,label=list(all_MC.keys()),bins=bins,range=(xmin,xmax),color=colors)
            signal_n, signal_bins = np.histogram(all_MC["ttHTobb"],weights=all_weight["weight_nominal"]["ttHTobb"],bins=bin_edges)
            ttdilepton_n, ttdilepton_bins = np.histogram(all_MC["$t\\bar{t}$ Dilepton"],weights=all_weight["weight_nominal"]["$t\\bar{t}$ Dilepton"],bins=bin_edges)
            ttnondilepton_n, ttnondilepton_bins = np.histogram(all_MC["$t\\bar{t}$ NonDilepton"],weights=all_weight["weight_nominal"]["$t\\bar{t}$ NonDilepton"],bins=bin_edges)
            signal_mcstaterr2 = np.histogram(all_MC["ttHTobb"],weights=all_weight["weight_nominal"]["ttHTobb"]**2,bins=bin_edges)[0]
            bg_mcstaterr2 = 0
            tt_n = ttdilepton_n + ttnondilepton_n
            bg_n = n[-1] - signal_n
            bg_nott = n[-1] - signal_n - tt_n
            signal_n = signal_n / ak.sum(signal_n)
            bg_n = bg_n / ak.sum(bg_n)
            bg_nott = bg_nott / ak.sum(bg_nott)
            tt_n = tt_n / ak.sum(tt_n)
            
            n_sys = {}
            for key in all_weight:
                if "Down" in key:
                    continue
                downString = key[:-2]+"Down"
                err2_up_total = np.zeros_like(n[-1])
                err2_down_total = np.zeros_like(n[-1])
                for label in all_weight[key]:
                    if key=="weight_nominal":
                        mcstat_err2 = np.histogram(all_MC[label], bins=bin_edges, weights=all_weight[key][label]**2)[0]
                        err2_up_total += mcstat_err2
                        err2_down_total += mcstat_err2
                        bg_mcstaterr2 += mcstat_err2
                    else:
                        hist_n_up, hist_bins_up = np.histogram(all_MC[label],weights=all_weight[key][label],bins=bin_edges)
                        hist_n_down, hist_bins_down = np.histogram(all_MC[label],weights=all_weight[downString][label],bins=bin_edges)
                        hist_n, hist_bins = np.histogram(all_MC[label],weights=all_weight["weight_nominal"][label],bins=bin_edges)
                        err_up = hist_n_up - hist_n
                        err_down = hist_n_down - hist_n
                        # Compute the flags to check which of the two variations (up and down) are pushing the nominal value up and down
                        up_is_up = err_up > 0
                        down_is_down = err_down < 0
                        # Compute the flag to check if the uncertainty is one-sided, i.e. when both variations are up or down
                        is_onesided = up_is_up ^ down_is_down
                        # Sum in quadrature of the systematic uncertainties taking into account if the uncertainty is one- or double-sided
                        err2_up_twosided = np.where(up_is_up, err_up**2, err_down**2)
                        err2_down_twosided = np.where(up_is_up, err_down**2, err_up**2)
                        err2_max = np.maximum(err2_up_twosided, err2_down_twosided)
                        err2_up_onesided = np.where(is_onesided & up_is_up, err2_max, 0)
                        err2_down_onesided = np.where(is_onesided & down_is_down, err2_max, 0)
                        err2_up_combined = np.where(is_onesided, err2_up_onesided, err2_up_twosided)
                        err2_down_combined = np.where(
                            is_onesided, err2_down_onesided, err2_down_twosided
                        )
                        err2_up_total += err2_up_combined
                        err2_down_total += err2_down_combined
                n_sys[key] = err2_up_total
                n_sys[downString] = err2_down_total
            ax.tick_params(axis='x', labelsize=22)
            ax.tick_params(axis='y', labelsize=22)
            ax.set_xlim(bin_edges[0],bin_edges[-1])
            #systematic error band
            err2_up = np.zeros_like(n[-1])
            err2_down = np.zeros_like(n[-1])
            for key in n_sys:
                if "Down" in key:
                    continue
                downString = key[:-2]+"Down"
                err2_up += n_sys[key]
                err2_down += n_sys[downString]
            up = n[-1] + np.sqrt(err2_up)
            down = n[-1] - np.sqrt(err2_down)
            ratio_up = np.where(n[-1] != 0, up / n[-1], 1)
            ratio_down = np.where(n[-1] != 0, down / n[-1], 1)
            unc_band = np.array([ratio_down, ratio_up])
            rax.fill_between(
                bin_edges,
                np.r_[unc_band[0], unc_band[0, -1]],
                np.r_[unc_band[1], unc_band[1, -1]],
                label="syst. unc.",
                color= [0.,0.,0.,0.4],
                facecolor= [0.,0.,0.,0.],
                hatch= "////",
                linewidth= 0,
                step= 'post',
                zorder= 2
            )
            rax.hlines(
                1.0, *bin_edges[[0, -1]], colors='gray', linestyles='dashed'
            )

            ax.legend(ncol=2,loc="upper right",fontsize=16)
            plt.xlabel(xlabel,fontsize=16)
            nData = [0.0]*len(bin_edges[:-1])
            for sample in dataSamples:
                data = p[sample][cat][var]
                if(quantity!="N" and col!="events"): 
                    pt_data = p[sample][cat][col+"_pt"]
                    sortIndices = ak.argsort(pt_data,ascending=False)
                    data = data[sortIndices]
                    if "pt" in modifier:
                        index = int(modifier.split("_")[1]) - 1
                        data = ak.mask(data, ak.num(data) > index)[:, index]
                    elif modifier=="max":
                        data = ak.max(data,axis=1)
                none_mask = ak.is_none(data)
                data = data[~none_mask]
                inf_mask = np.isinf(data)
                data = data[~inf_mask]
                if(modifier=="all" and quantity!="N" and col!="events"):
                    data = ak.flatten(data)
                n_data, bin_edges_data, patches_data = a2.hist(data,bins=bin_edges)
                nData += n_data
            error = np.sqrt(nData)
            binCenters = (bin_edges_data[:-1] + bin_edges_data[1:]) / 2
            ratio = nData/n[-1]
            ratio_error = error/n[-1]
            ratio_error[np.isnan(ratio_error)] = np.inf
            ax.errorbar(binCenters, nData, yerr=error,linestyle='solid',markersize=5.0,color='black',elinewidth=1,linewidth=0,marker='.')
            #ax.set_yscale("log")
            rax.errorbar(binCenters, ratio, yerr=ratio_error,linestyle='solid',markersize=5.0,color='black',elinewidth=1,linewidth=0,marker='.')
            rax.set_ylim((0.5, 1.5))
            rax.yaxis.set_label_coords(-0.075, 1)
            rax.tick_params(axis='x', labelsize=22)
            rax.tick_params(axis='y', labelsize=22)
            max_data = max(nData)
            max_MC = max(n[-1])
            if(max_data>max_MC):
                ax.set_ylim((0, 2.0 * max_data))
            else:
                ax.set_ylim((0, 2.0 * max_MC))
        else:
            f1, a1 = plt.subplots()
            if(bins==-1):
                n, bin_edges, patches = a1.hist(list(all_MC.values()),weights=list(all_weight["weight_nominal"].values()),stacked=True,label=list(all_MC.keys()))
            else:
                n, bin_edges, patches = a1.hist(list(all_MC.values()),weights=list(all_weight["weight_nominal"].values()),stacked=True,label=list(all_MC.keys()),bins=bins,range=(xmin,xmax))
            ax.legend(ncol=2,loc="upper right")

        plt.close(f2)
        if "pt" in modifier:
            stringMod = modifier.split("_")[1]
        else:
            stringMod = modifier
        filepath = f"hists_optimizedLepIsoandbbTagging_0GeVMET/2018/{cat}/{col}"
        if not os.path.exists(filepath):
            os.makedirs(filepath)
        plt.savefig(filepath+"/"+var+"_"+stringMod+".pdf")
        plt.show()
        plt.clf()
        plt.cla()
        #normalized signal vs bg
        binCenters = (bin_edges[:-1] + bin_edges[1:]) / 2
        bg_mcstaterr2 = bg_mcstaterr2 - signal_mcstaterr2
        plt.figure(1)
        ax1 = plt.subplot(3,1,(1,2))
        signal_n = np.append(signal_n,signal_n[-1])
        bg_n = np.append(bg_n,bg_n[-1])
        plt.step(bin_edges,signal_n,color="red",label="signal",where="post")
        plt.step(bin_edges,bg_n,color="blue",label="bkg",where="post")
        ax1.legend(ncol=2,loc="upper right",fontsize=16)
        ax2 = plt.subplot(3,1,3,sharex=ax1)
        plt.scatter(binCenters,signal_n[:-1]/bg_n[:-1],color="black",marker=".")
        ax2.hlines(
                1.0, *bin_edges[[0, -1]], colors='gray', linestyles='dashed'
            )
        ax2.set_yscale("log")
        plt.xlabel(xlabel,fontsize=16)
        plt.show()
        plt.close()
        #normalized signal vs tt vs other bg
        plt.figure(2)
        ax1 = plt.subplot(3,1,(1,2))
        tt_n = np.append(tt_n,tt_n[-1])
        bg_nott = np.append(bg_nott,bg_nott[-1])
        plt.step(bin_edges,signal_n,color="red",label="signal",where="post")
        plt.step(bin_edges,bg_nott,color="blue",label="other bkg",where="post")
        plt.step(bin_edges,tt_n,color="green",label="tt bkg",where="post")
        ax1.legend(ncol=2,loc="upper right",fontsize=16)
        ax2 = plt.subplot(3,1,3,sharex=ax1)
        plt.scatter(binCenters,signal_n[:-1]/bg_n[:-1],color="black",marker=".")
        ax2.hlines(
                1.0, *bin_edges[[0, -1]], colors='gray', linestyles='dashed'
            )
        ax2.set_yscale("log")
        plt.xlabel(xlabel,fontsize=16)
        plt.show()
        plt.close()

'''

# muonIso: miniPFRelIso_chg WP: 0.4 eleIso: miniPFRelIso_chg WP: 0.4 bbWP: 0.0399
ElectronIsoWP = [0.4]
MuonMiniIsoWP = [0.4]
MuonPFIsoWP = [0.4]
BBWorkingPoints = [0.0399]
ElectronIsoCollections = ["miniPFRelIso_chg" ]
MuonIsoCollections =  ["miniPFRelIso_chg"]

metWPs = [0, 5, 10, 15, 20, 25, 30, 35, 40, 50] # 0-40 every 5 GeV # dilepton mass cut, custom cut function.py commented out,
# below are whatwe had originally
#ElectronIsoWP = [0.06]
#MuonPFIsoWP = [0.25]
#BBWorkingPoints = [0.0]
#ElectronIsoCollections = ["pfRelIso03_all" ]
#MuonIsoCollections =  ["pfRelIso04_all"]
#Note, the newest paper from BTV on bbtagging references WPs but doesn't seem to list them
#ElectronIsoCollections = ["jetRelIso", "pfRelIso03_all","miniPFRelIso_all", "pfRelIso03_chg","miniPFRelIso_chg", "dr03EcalRecHitSumEt", "dr03HcalDepth1TowerSumEt", "dr03TkSumPt", "dr03TkSumPtHEEP"]
#MuonIsoCollections = ["jetRelIso","pfRelIso03_all", "miniPFRelIso_all" , "pfRelIso03_chg","miniPFRelIso_chg", "pfRelIso04_all"]
f = open("METPt_Optimization_0-50_usingisofromRyansStudy_April4.txt", "a")

#metWP = 0.
jetWPs = [200., 225., 250., 275., 300., 325., 350.]
jetWP = 200.
for metWP in metWPs:
#for jetWP in jetWPs:
    for muonIso in MuonIsoCollections:
        if "mini" in muonIso:
            muonWPs = MuonMiniIsoWP
        else:
            muonWPs = MuonPFIsoWP
        for eleIso in ElectronIsoCollections:
            for muonWP in muonWPs:
                for eleWP in ElectronIsoWP:
                    for bbWP in BBWorkingPoints:
                        p = {}
                        parquetDir = "/afs/cern.ch/work/a/asparker/public/ttHboosted_March20_merge"
                        cutFlow = {"Signal" : {"baseline": 0, "ee": 0, "emu": 0, "mumu": 0}, "TT" : {"baseline": 0, "ee": 0, "emu": 0, "mumu": 0}, "OtherBkg" : {"baseline": 0, "ee": 0, "emu": 0, "mumu": 0}}
                        for sample in os.listdir(parquetDir):
                            if(("TTbb" in sample) or ("DATA" in sample)):
                                continue
                            #print(sample,"muonIso",muonIso,"WP",muonWP,"eleIso",eleIso,"WP",eleWP,"bbWP",bbWP)
                            p[sample] = {}
                            #for cat in os.listdir(parquetDir+"/"+sample):
                            for cat in ["baseline"]:
                                p[sample][cat] = ak.from_parquet(parquetDir+"/"+sample+"/"+cat+"/nominal/merge.parquet")
                                eleMask = p[sample][cat]["ElectronGood_"+eleIso]<eleWP
                                muonMask = p[sample][cat]["MuonGood_"+muonIso]<muonWP
                                metMask = p[sample][cat]["MET_pt"]>metWP
                               
                                if("v7" in sample):
                                    fatMask = p[sample][cat]["FatJetGood_btagDDBvL"]>bbWP
                                else:
                                    fatMask = p[sample][cat]["FatJetGood_btagDDBvLV2"]>bbWP

                                #metPt  = p[sample][cat]["MET_pt"][metMask]
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
                                #metPt = metPt[mask_lepton_cleaning] 
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
                                eventMask = (nlep==2) & (ak.sum(lepPt>=25.0,axis=1)>=2) & (njet>=2) & (ll_charge==0) & (mass_cut) & (ak.sum(fatPt>=200,axis=1)>=1) & metMask  #(metPt[0] >= metWP)  #(ak.sum(metPt>=metWP,axis=0))
                                #print("metPt[0] ", metPt[0])
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
                            print("muonIso",muonIso,"WP",muonWP,"eleIso",eleIso,"WP",eleWP,"bbWP",bbWP, "metWP", metWP, "jetWP", jetWP )
                            f.write("muonIso: {} WP: {} eleIso: {} WP: {} bbWP: {} metWP: {} jetWP: {}\n".format(muonIso,muonWP,eleIso,eleWP,bbWP,metWP,jetWP))
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








'''
