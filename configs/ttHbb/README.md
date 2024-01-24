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
```
# Install in EDITABLE mode
```
pip install -e .[dev]
```

### Clone locally the AnalysisConfigs repo
```
git clone git@github.com:ryanm124/AnalysisConfigs.git
```

## To run the code

### make datasets
```
build_datasets.py --cfg datasets/MC_all.json -o
```
#### run config with dataset names (configs 1-9 and 12-20 are the ones you need to add up to get all data and MC for 2018)
```
runner.py --cfg newconfig1.py  -o out_Jan24
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
