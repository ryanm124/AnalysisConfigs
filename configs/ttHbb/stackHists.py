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
import math



filename = "/afs/cern.ch/work/r/rmccarth/public/output_allNoDupes.coffea"


o = load(filename)
'''
print("o.keys()")
print(o.keys())
print("o['variables'].keys()")
print(o['variables'].keys())
print("o['columns'].keys()")
print(o['columns'].keys())

print("o['datasets_metadata'].keys()")
print(o['datasets_metadata'].keys())
'''



print("CUTFLOW")
print(o['cutflow'].keys())

'''
print("o['sumw']['baseline']")
print(o['sumw']['baseline'])

print("o['sumw']['baseline']['TTToSemiLeptonic__2018']")
print(o['sumw']['baseline']['TTToSemiLeptonic__2018'])
'''

# selection                                                                                                                                    
selection = 'mumu'  #'emu'#'ee'#'baseline'

vars = {}
for i, var_name in enumerate(o['variables'].keys()):
    hists = {}
    vars[var_name] = hists
    meras = 0
    eeras = 0
    for sam_name in o['columns'].keys():
        sam_namey = sam_name[:-9] + '__2018'
        print("sam_namey")
        print(sam_namey)
        if 'ttHToNonbb' in sam_name: sam_namey ='ttHToNonbb_Powheg_2018'
        if 'DYJetsToLL_M-50' in sam_name: sam_namey ='DYJetsToLL_M-50_v7__2018'
        if 'ZJetsToQQ' in sam_name: sam_namey = sam_name[:-9] + '_v7__2018'
        if 'ST' in sam_name: sam_namey = sam_name[:-9] + '_v7__2018'
        if 'QCD' in sam_name: sam_namey = sam_name[:-9] + '_v7__2018'
        if 'WJets' in sam_name: sam_namey = sam_name[:-9] + '_v7__2018' 
        if 'WJetsToLNu' in sam_name: sam_namey = sam_name[:-9] + '__2018'
        if 'TTWJets' in sam_name: sam_namey = sam_name[:-9] + '__2018'
        if 'TTGJets' in sam_name: sam_namey = sam_name[:-9] + '__2018'
        if 'THW' in sam_name: sam_namey = sam_name[:-9] + '__2018'
        if 'WW' in sam_name: sam_namey = sam_name[:-9] + '__2018'
        if 'WZ' in sam_name: sam_namey = sam_name[:-9] + '__2018'
        if 'ZZ' in sam_name: sam_namey = sam_name[:-9] + '__2018'
        if 'ttHTobb' in sam_name: sam_namey = sam_name[:-9] + '_2018'
        print("NEW sam_name")
        print(sam_name)
        print(var_name)
        if 'SingleMuon' in sam_name: 
            if meras == 0:
                sam_namey = sam_name[:-9] + '_2018_EraA'
            if meras == 1:
                sam_namey = sam_name[:-9] + '_2018_EraB'
            if meras == 2:
                sam_namey = sam_name[:-9] + '_2018_EraC'
            if meras == 3:
                sam_namey = sam_name[:-9] + '_2018_EraD'

            meras+=1
        if 'EGamma' in sam_name:
            if eeras == 0:
                sam_namey = sam_name[:-9] + '_2018_EraA'
            if eeras == 1:
                sam_namey = sam_name[:-9] + '_2018_EraB'
            if eeras == 2:
                sam_namey = sam_name[:-9] + '_2018_EraC'
            if eeras == 3:
                sam_namey = sam_name[:-9] + '_2018_EraD'

            eeras+=1

        varHist = o['variables'][var_name][sam_name[:-9]][sam_namey ]

        # baseline is 0th hist
        
        aindex = 0
        if selection == 'baseline': aindex = 0
        if selection ==	'ee': aindex = 1
        if selection == 'emu': aindex = 2
        if selection == 'mumu': aindex = 3
        h = varHist.stack("cat").project(varHist.axes[-1].name)[aindex]
        hists[sam_name[:-9]]=h

import yaml

# Load the data from the YAML file
with open('/afs/cern.ch/user/a/asparker/public/dec_ttH/PocketCoffea/AnalysisConfigs/configs/ttHbb/params/plotting_style.yaml', 'r') as yaml_file:
    plotting_styles = yaml.safe_load(yaml_file)


# selection 
#selection = 'baseline'

# Assuming your cutflow data is stored in the 'cutflow_data' dictionary
cutflow_data = o['cutflow'][selection]  #['mumu']#['baseline']
uncerts = o['sumw2'][selection] # ['mumu']



print(o['sum_genweights'].keys())
genweights = o['sum_genweights']

cutflow_data = o['cutflow'][selection] 

columns = o['columns']

# Define LaTeX table header
latex_table = "\\begin{table}[htbp]\n"
latex_table += "\\centering\n"
latex_table += "\\begin{tabular}{|l|c|}\n"
latex_table += "\\hline\n"
latex_table += "Process & Count \\\\\n"
latex_table += "\\hline\n"

# Iterate through the groups in the plotting_styles YAML and add rows to the LaTeX table
for group, processes in plotting_styles['plotting_style']['samples_groups'].items():
    group_count = 0
    group_uncert = 0
    gCount = []
    gUnc = [] 
    print("processes")
    print(processes)
    for process in processes:
        print("process")
        print(process)
        process2 = process + "__2018"
        if 'ttHToNonbb' in process2: process2 ='ttHToNonbb_Powheg_2018'
        if 'DYJetsToLL_M-50' in process2: process2 ='DYJetsToLL_M-50_v7__2018'
        if 'ZJetsToQQ' in process2: process2 = process+ '_v7__2018'
        if 'ZJetsToQQ' in process2: process2 = process+ '_v7__2018'
        if 'ST' in process2: process2 = process+ '_v7__2018'
        if 'QCD' in process2: process2 = process+ '_v7__2018'
        if 'ttHTobb' in process2: process2 = process+ '_2018'
        if 'WJetsToQQ' in process2 : process2 = process+ '_v7__2018'
        if 'TTWJetsToQQ'  in process2: process2 = process+ '__2018'
        
        cutflow_weight = 0
 
        weights1 = ak.from_numpy(columns[process+'__nominal'][process2][selection]['weight_nominal'].value)
        weights = ak.sum(weights1,axis=0)

        squared_errors_sum = 0
        sum_genweights = genweights[process2]
        for w in zip( weights1 ):
            # Calculate weight / sum_genweight for each event
            ratio = w / sum_genweights

            # Square the result and add to the sum
            squared_errors_sum += ratio ** 2

        # Take the square root of the sum
        error2 = np.sqrt(squared_errors_sum)

        
        
        #print("Resulting error2:", error2)
 
        cutflow_weight = weights
        print("cutflow_weight")
        print(cutflow_weight)
        print("sum_genweights")
        print(sum_genweights)

        if process2 in cutflow_data:
            print(process)
            gCount.append( cutflow_weight / sum_genweights    )
            gUnc.append( math.sqrt( uncerts[process2][process] ) )
    print(gCount)
    print(gUnc)
    for c in gCount: 
        group_count += c
    for u in gUnc:
        group_uncert += u**2
        
    # Replace group names with the ones from the YAML labels_mc section
    group_label = plotting_styles['plotting_style']['labels_mc'].get(group, group)
    
    group_uncert2 = math.sqrt(group_uncert)

    latex_table += f"{group_label} & {group_count:.0f} \pm {group_uncert2:.2f}  \\\\\n"

# Add LaTeX table footer
latex_table += "\\hline\n"
latex_table += "\\end{tabular}\n"
latex_table += "\\caption{Baseline Cutflow table}\n"
latex_table += "\\label{tab:cutflow}\n"
latex_table += "\\end{table}"

# Print the LaTeX table
print(latex_table)
