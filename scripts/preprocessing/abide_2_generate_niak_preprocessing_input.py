import numpy as np
import pandas as pd
import pathlib as pal
from scipy import io as sio


root_p = pal.Path('/project/6003287/DATA/ABIDE_2/')
raw_p = root_p / 'RAW'
pheno_p = root_p / 'Pheno/ABIDEII_Composite_Phenotypic.csv'

pheno = pd.read_csv(pheno_p, encoding="ISO-8859-1")
# Rename ETH to ETHZ
pheno.loc[pheno['SITE_ID'] == "ABIDEII-ETH_1", ['SITE_ID']] = 'ABIDEII-ETHZ_1'
# Rename OILH to ONRC
pheno.loc[pheno['SITE_ID'] == "ABIDEII-OILH_2", ['SITE_ID']] = 'ABIDEII-ONRC_2'
exclude_these_subjects = ['sub0029007', 'sub0029008', 'sub0029009', 'sub0029010', 'sub0028947',
                          'sub0028954', 'sub0030172', 'sub0028685', 'sub0028697', 'sub0028727']

for site in pheno.SITE_ID.unique():
    site_name = site.replace('-', '_')
    files_in = dict()
    site_p = raw_p / site
    subjects = pheno.query(f'SITE_ID == "{site}"')['SUB_ID'].unique()
    for sub in subjects:
        sub_name = f'sub{sub:07}'
        if sub_name in exclude_these_subjects:
            print(f'I found {sub_name} from {site} in the list of to be excluded subjects. Will exclude.')
            continue
        sub_p = site_p / str(sub)
        ses_p = sub_p / 'session_1'
        anat = list((ses_p / 'anat_1').glob('*.nii.gz'))
        rest = list(ses_p.glob('rest_*/*.nii.gz'))
        if not anat or not rest:
            continue
        files_in[sub_name] = dict()
        files_in[sub_name]['anat'] = str(anat[0])
        files_in[sub_name]['fmri'] = dict()
        files_in[sub_name]['fmri']['session1'] = dict()
        for rest_p in rest:
            run_name = rest_p.parent.name.replace('_', '')
            files_in[sub_name]['fmri']['session1'][run_name] = str(rest_p)

    # Build the options
    opt = dict()
    opt['folder_out'] = str(root_p / 'PREPROCESS_NIAK' / site)  # Path to preprocessed data for this site
    opt['size_output'] = 'quality_control'

    opt['slice_timing'] = dict()
    opt['slice_timing']['type_acquisition'] = 'interleaved ascending'  # slice acquisition order
    opt['slice_timing']['type_scanner'] = 'Siemens'  # Manufacturer
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
    opt['psom']['max_queued'] = int(np.floor(len(subjects) * 0.2)) + 1

    variables = dict()
    variables['files_in'] = files_in
    variables['opt'] = opt
    sio.savemat(str(root_p / 'Support_NIAK' / f'{site_name}_niak.mat'), variables)
