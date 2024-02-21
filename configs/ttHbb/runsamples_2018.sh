#!/bin/bash

# Define an array to store the unique sample names
declare -a samples=(
    "ttHTobb"
    "ttHToNonbb"
    "TTbbDiLeptonic_Powheg"
    "TTbbSemiLeptonic_Powheg"
    "TTbbHadronic_Powheg"
    "THW"
    "TTToSemiLeptonic_2018"
    "TTTo2L2Nu"
    "TTToHadronic"
    "TTZToQQ"
    "TTZtoLLNuNu"
    "TTWJetsToQQ"
    "TTWJetsToLNu"
    "TTGJets"
    "WW"
    "ZZ"
    "WZ"
    "WJetsToLNu_HT-200To400"
    "WJetsToLNu_HT-400To600"
    "WJetsToLNu_HT-600To800"
    "WJetsToLNu_HT-800To1200"
    "WJetsToLNu_HT-1200To2500"
    "WJetsToLNu_HT-2500ToInf"
    "WJetsToQQ_HT400to600_v7"
    "WJetsToQQ_HT600to800_v7"
    "WJetsToQQ_HT800toInf_v7"
    "ZJetsToQQ_HT400to600_v7"
    "ZJetsToQQ_HT600to800_v7"
    "ZJetsToQQ_HT800toInf_v7"
    "DYJetsToLL_M-50_v7"
    "ST_s-channel_4f_leptonDecays_v7"
    "ST_t-channel_top_4f_InclusiveDecays_v7"
    "ST_t-channel_antitop_4f_InclusiveDecays_v7"
    "ST_tW_top_5f_inclusiveDecays_v7"
    "ST_tW_antitop_5f_inclusiveDecays_v7"
    "QCD_HT500to700_v7"
    "QCD_HT700to1000_v7"
    "QCD_HT1000to1500_v7"
    "QCD_HT1500to2000_v7"
    "QCD_HT2000toInf_v7"
)

# Loop through unique samples and execute runner.py command
for sample in "${samples[@]}"; do
    echo "Running command for sample: $sample"
    runner.py --cfg config_general.py -y 2018 -sa "$sample" -o production_Feb21 --full
done
