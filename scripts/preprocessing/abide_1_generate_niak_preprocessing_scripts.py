import json
import pandas as pd
import pathlib as pal

root_p = pal.Path('/home/surchs/sim_big/DATA/ABIDE_1/')
script_paths = pal.Path('/home/surchs/local_projects/abide_univariate/scripts/preprocessing/abide_1/')
raw_p = root_p / 'RAW'

script_str = '''
%%% NKI TRT preprocessing pipeline
% Script to run a preprocessing pipeline analysis on the HCP database.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting input/output files %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
root_p = '/project/6003287/DATA/ABIDE_1';
site_name = "{}";
load([root_p filesep 'Support_NIAK' filesep site_name '_niak.mat']);
opt.flag_test = false;
[pipeline,opt] = niak_pipeline_fmri_preprocess(files_in,opt);
'''

sites = [p for p in raw_p.glob('*') if p.is_dir()]
sites_json = [site / 'task-rest_bold.json' for site in sites]

# Find the scanning parameters
data = {key: list() for key in ['Site', 'Manufacturer', 'SliceAcquisitionOrder', 'TR', 'Path']}
for site_id, jp in enumerate(sites_json):
    site_name = jp.parent.name
    data['Site'].append(site_name)
    data['Path'].append(sites[site_id])
    with open(jp, 'r') as f:
        text = json.load(f)
    if 'Manufacturer' in text.keys():
        data['Manufacturer'].append(text['Manufacturer'])
    else:
        print(site_name, text.keys())
        data['Manufacturer'].append('unknown')

    if 'SliceAcquisitionOrder' in text.keys():
        data['SliceAcquisitionOrder'].append(text['SliceAcquisitionOrder'])
    else:
        print(site_name, text.keys())
        data['SliceAcquisitionOrder'].append('unknown')
    if 'RepetitionTime' in text.keys():
        data['TR'].append(text['RepetitionTime'])
    else:
        print(site_name, text.keys())
        data['TR'].append('unknown')
table = pd.DataFrame(data)

for site_id, site in table.iterrows():
    # See if we have a site without SliceTimeOrder
    if site['SliceAcquisitionOrder'] == 'unknown':
        print(f"Site {site['Site']} has an unknown acquisition order. I am still including this because we don't "
              f"do slice time correction.")
    site_name = site['Site']
    with open(script_paths / f'{site_name}_niak_preproc.m', 'w') as f:
        f.write(script_str.format(site_name))
