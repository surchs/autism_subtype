import json
import pandas as pd
import pathlib as pal


root_p = pal.Path('/home/surchs/data/surchs/ABIDE/ABIDE_2/')
script_paths = pal.Path('/home/surchs/Projects/abide_univariate/scripts/preprocessing/abide_2/')
raw_p = root_p / 'RAW'

script_str = '''
%%% NKI TRT preprocessing pipeline
% Script to run a preprocessing pipeline analysis on the HCP database.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting input/output files %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
root_p = '/project/6003287/DATA/ABIDE_2';
site_name = "{}";
load([root_p filesep 'Support_NIAK' filesep site_name '_niak.mat']);
opt.flag_test = false;
[pipeline,opt] = niak_pipeline_fmri_preprocess(files_in,opt);
'''

root_p  = pal.Path('/home/surchs/data/surchs/ABIDE/ABIDE_2')
pheno_p = root_p / 'Pheno/ABIDEII_Composite_Phenotypic.csv'

pheno = pd.read_csv(pheno_p, encoding="ISO-8859-1")
# Rename ETH to ETHZ
pheno.loc[pheno['SITE_ID'] == "ABIDEII-ETH_1", ['SITE_ID']] = 'ABIDEII-ETHZ_1'
# Rename OILH to ONRC
pheno.loc[pheno['SITE_ID'] == "ABIDEII-OILH_2", ['SITE_ID']] = 'ABIDEII-ONRC_2'


for site in pheno.SITE_ID.unique():
    site_name = site.replace('-', '_')
    with open(script_paths / f'{site_name}_niak_preproc.m', 'w') as f:
        f.write(script_str.format(site_name))
