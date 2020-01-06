import re
import json
import numpy as np
import pandas as pd
import pathlib as pal
from scipy import io as sio


root_p = pal.Path('/project/6003287/DATA/ABIDE_1/')
raw_p = root_p / 'RAW'
sites = [p for p in raw_p.glob('*') if p.is_dir()]
sites_json = [site / 'task-rest_bold.json' for site in sites]
exclude_these_subjects = ['sub0050566']
um1_exclusions = ['sub0050294', 'sub0050295', 'sub0050298', 'sub0050299',
                  'sub0050313', 'sub0050326', 'sub0050333', 'sub0050335', 'sub0050336', 'sub0050337', 'sub0050338',
                  'sub0050344', 'sub0050346', 'sub0050350', 'sub0050363', 'sub0050368', 'sub0050372', 'sub0050377']
um2_exclusions = ['sub0050402', 'sub0050403',
                  'sub0050404', 'sub0050405', 'sub0050407', 'sub0050408',
                  'sub0050411', 'sub0050413', 'sub0050414', 'sub0050415',
                  'sub0050419', 'sub0050423', 'sub0050425', 'sub0050427', 'sub0050421']
exclude_these_subjects += um1_exclusions
exclude_these_subjects += um2_exclusions

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
    site_name = site['Site']
    if site['SliceAcquisitionOrder'] == 'unknown':
        print(f"Site {site['Site']} has an unknown acquisition order. I am still including this because we don't "
              f"do slice time correction.")
    # Find subjects
    site_path = site['Path']
    subjects = list(site_path.glob('sub*'))

    # Find all subjects
    n_runs = 0
    files_in = dict()
    for subject in subjects:
        subject_name = subject.name.replace('-', '')
        if subject_name in exclude_these_subjects:
            print(f'I found {subject_name} from {site_name} in the list of to be excluded subjects. Will exclude.')
            continue
        anat_hits = list(subject.glob('anat/*T1w.nii.gz'))
        if len(anat_hits) > 1:
            print(f'{subject_name} has more than one T1: {anat_hits}')
            continue
        files_in[subject_name] = dict()
        files_in[subject_name]['anat'] = str(anat_hits[0])
        # Find runs
        runs = list(subject.glob('func/*bold.nii.gz'))
        files_in[subject_name]['fmri'] = dict()
        files_in[subject_name]['fmri']['session1'] = dict()
        for run in runs:
            run_name = re.search(r'run-\d', run.name).group().replace('-', '')
            files_in[subject_name]['fmri']['session1'][run_name] = str(run)
            n_runs += 1

    # Build the options
    opt = dict()
    opt['folder_out'] = str(root_p / 'PREPROCESS_NIAK' / site['Site'])  # Path to preprocessed data for this site
    opt['size_output'] = 'quality_control'

    opt['slice_timing'] = dict()
    opt['slice_timing']['type_acquisition'] = site['SliceAcquisitionOrder'].lower()  # slice acquisition order
    opt['slice_timing']['type_scanner'] = site['Manufacturer']  # Manufacturer
    opt['slice_timing']['delay_in_tr'] = 0.0
    opt['slice_timing']['suppress_vol'] = 4
    opt['slice_timing']['flag_nu_correct'] = 0
    opt['slice_timing']['arg_nu_correct'] = '-distance 200'
    opt['slice_timing']['flag_center'] = 0
    opt['slice_timing']['flag_skip'] = 1

    opt['motion'] = dict()
    opt['motion']['session_ref'] = 'session1'

    opt['resample_vol'] = dict()
    opt['resample_vol']['interpolation'] = 'trilinear'
    opt['resample_vol']['voxel_size'] = [3, 3, 3]
    opt['resample_vol']['flag_skip'] = 0

    opt['t1_preprocess'] = dict()
    opt['t1_preprocess']['nu_correct'] = dict()
    opt['t1_preprocess']['nu_correct']['arg'] = '-distance 75'

    opt['time_filter'] = dict()
    opt['time_filter']['hp'] = 0.01
    opt['time_filter']['lp'] = np.inf

    opt['regress_confounds'] = dict()
    opt['regress_confounds']['flag_slow'] = True
    opt['regress_confounds']['flag_high'] = False
    opt['regress_confounds']['flag_motion_params'] = True
    opt['regress_confounds']['flag_pca_motion'] = True
    opt['regress_confounds']['pct_var_explained'] = 0.95

    opt['regress_confounds'][
        'flag_wm'] = True  # Turn on/off the regression of the average white matter signal (true: apply / false : don't apply)
    opt['regress_confounds'][
        'flag_vent'] = True  # Turn on/off the regression of the average of the ventricles (true: apply / false : don't apply)
    opt['regress_confounds'][
        'flag_gsc'] = False  # Turn on/off the regression of the PCA-based estimation of the global signal (true: apply / false : don't apply)
    opt['regress_confounds'][
        'flag_scrubbing'] = True  # Turn on/off the scrubbing of time frames with excessive motion (true: apply / false : don't apply)
    opt['regress_confounds'][
        'thre_fd'] = 0.4  # The threshold on frame displacement that is used to determine frames with excessive motion in the

    opt['smooth_vol'] = dict()
    opt['smooth_vol']['fwhm'] = 6  # Full-width at maximum (FWHM) of the Gaussian blurring kernel, in mm.
    opt['smooth_vol']['flag_skip'] = 0  # Skip spatial smoothing (0: don't skip, 1 : skip)

    opt['psom'] = dict()
    opt['psom']['qsub_options'] = '--account rpp-aevans-ab --time=00-05:00 --ntasks=1 --cpus-per-task=3 --mem-per-cpu=4G'
    opt['psom']['max_queued'] = int(np.floor(n_runs * 0.2))


    variables = dict()
    variables['files_in'] = files_in
    variables['opt'] = opt
    sio.savemat(str(root_p / 'Support_NIAK' / f'{site_name}_niak.mat'), variables)
