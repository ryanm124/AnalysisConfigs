# ttH Boosted 

***NOTE : This analysis uses the [NANOAOD](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD) data format. ***


## Set Up the Code

You need a grid proxy before entering the singularity container
```
voms-proxy-init -voms cms -rfc --valid 168:0

```
### Enter the singularity

```
apptainer shell -B ${XDG_RUNTIME_DIR} -B /afs -B /cvmfs/cms.cern.ch --bind /tmp  -B /etc/grid-security/ -B /etc/vomses --bind /eos/cms/  --bind /etc/sysconfig/ngbauth-submit  --env KRB5CCNAME=${XDG_RUNTIME_DIR}/krb5cc /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-cc7-latest
```


### Create a local virtual environment using the packages defined in the singularity image
```
python -m venv --system-site-packages myenv
```

### Activate the environment
```
source myenv/bin/activate
```
### Clone locally the PocketCoffea repo
```
git clone git@github.com:ryanm124/PocketCoffea.git
cd PocketCoffea
git checkout SRstudies
```
# Install in EDITABLE mode
```
pip install -e .[dev]
```

### Clone locally the AnalysisConfigs repo
```
git clone git@github.com:ryanm124/AnalysisConfigs.git
cd AnalysisConfigs
git checkout SRstudies
```

## To run the code

### make datasets
```
build_datasets.py --cfg datasets/MC_all.json -bs T1_US_FNAL_Disk -o
build_datasets.py --cfg datasets/datasets_EGamma_definitions.json -bs  T1_US_FNAL_Disk  -o 
build_datasets.py --cfg datasets/datasets_SingleMuon_definitions.json -bs  T1_US_FNAL_Disk  -o

build_datasets.py --cfg datasets/QCD_500_2018.json -ws T1_DE_KIT_Disk -o 
build_datasets.py --cfg datasets/QCD_700_2018.json -ws T1_FR_CCIN2P3_Disk -o 

build_datasets.py --cfg datasets/TT_2018.json -ws T2_BE_IIHE -o 
build_datasets.py --cfg datasets/QCD_1000_2018.json -ws T2_AT_Vienna -o
build_datasets.py --cfg datasets/QCD_1500_2018.json -ws T1_DE_KIT_Disk -o

build_datasets.py --cfg datasets/QCD_1000_2018.json -ws T2_AT_Vienna -o

build_datasets.py --cfg datasets/WJets_2500_2018.json -ws T2_CH_CERN -o

```
#### run config with dataset names (configs 1-9 and 12-20 are the ones you need to add up to get all data and MC for 2018)
be sure to change the output file location for the parquet files and to creat that directory e.g. /eos/cms/store/user/asparker/ttHboosted_apri24/chunks
Note that you will need to move the files from EOS after processing them as the EOS files are not viewable in the jupyter notebooks we use within the signularity container
```
runner.py --cfg config_EGamma_data.py -o production_March14 --executor dask@lxplus --custom-run-options custom_run_options.yaml
runner.py --cfg config_SingleMuon_data_2018.py -o production_March14 --executor dask@lxplus --custom-run-options custom_run_options.yaml

runner.py --cfg config_QCD_2018.py -o production_March14 --executor dask@lxplus --custom-run-options custom_run_options.yaml
runner.py --cfg config_TTTo2L2Nu_2018.py -o production_March14 --executor dask@lxplus --custom-run-options custom_run_options.yaml
runner.py --cfg config_ST_2018.py -o production_March14 --executor dask@lxplus --custom-run-options custom_run_options.yaml
runner.py --cfg config_Tops_2018.py -o production_March14 --executor dask@lxplus --custom-run-options custom_run_options.yaml
runner.py --cfg config_ttH_2018.py -o production_March14 --executor dask@lxplus --custom-run-options custom_run_options.yaml
runner.py --cfg config_TTbb_Powheg_2018.py -o production_March14  --executor dask@lxplus --custom-run-options custom_run_options.yaml
runner.py --cfg config_bosons_2018.py -o production_March14  --executor dask@lxplus --custom-run-options custom_run_options.yaml


```
### merge outfiles
```
merge_outputs.py -i  output_newconfig7_all.coffea output_newconfig8_all.coffea output_newconfig3_all.coffea output_newconfig6_all.coffea    -o output_most_2018.coffea
```
### make plots with default PocketCoffea plotting script
```
make_plots.py --cfg configurator.pkl -o plots_most -i output_most_2018.coffea -j 6 --log -e jet_eta_pt_leading jet_eta_pt_all subleading_jetpt_MET  --overwrite
```
### make plots with Ryan's custom plotting script
```
jupyter HomebrewPlotting.ipynb
```
