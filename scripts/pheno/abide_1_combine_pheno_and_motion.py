import re
import pandas as pd
import pathlib as pal

raw_p = pal.Path('/home/surchs/data/surchs/ABIDE/ABIDE_1/')
root_p = pal.Path('../../data/')
prep_d = raw_p / 'PREPROCESS_NIAK/'
report_d = raw_p / 'Report'
pheno_p = root_p / 'pheno/ABIDE1_Pheno_prepared.tsv'
pheno_out_p = root_p / 'pheno/ABIDE1_Pheno_prepared_motion.tsv'

pheno = pd.read_csv(pheno_p, sep='\t')
# Some sites are removed because they didn't complete successfully
drop_sites = ['\"UM_1\"', '\"UM_2\"', '\"CMU_a\"', '\"CMU_b\"']
query = f'site_name != {" and site_name != ".join(drop_sites)}'
pheno_select = pheno.query(query)
site_select = list(pheno_select.site_name.unique())

total_list = list()
for site in site_select:
    if site == 'Stanford':
        continue
    site_report_p = report_d / site
    site_sub_p = site_report_p / 'assets/group/js/listSubject.js'
    site_run_p = site_report_p / 'assets/group/js/listRun.js'

    qc_p = prep_d / site / 'quality_control/group_motion/qc_scrubbing_group.csv'
    qc = pd.read_csv(qc_p)
    qc.rename(columns={' ': 'file_name'}, inplace=True)
    # Strip the white space
    qc.rename(columns={col: col.strip(' ') for col in qc.columns}, inplace=True)

    for rid, row in qc.iterrows():
        sub, ses, run = row.file_name.split('_')
        sid = int(re.search(r"(?<=sub00)\d+", sub).group())
        sub_dict = pheno.query(f'SUB_ID=={sid}').iloc[0].to_dict()
        sub_dict['session'] = int(re.search(r'(?<=session)\d', ses).group())
        sub_dict['run'] = int(re.search(r'(?<=run)\d', run).group())
        sub_dict['fd'] = row.FD
        sub_dict['fd_scrubbed'] = row.FD_scrubbed
        sub_dict['frames_ok'] = row.frames_OK
        sub_dict['frames_scrubbed'] = row.frames_scrubbed
        total_list.append(sub_dict)

columns = list(pheno.columns) + ['session', 'run', 'fd', 'fd_scrubbed', 'frames_ok', 'frames_scrubbed']
table = pd.DataFrame(total_list, columns=columns)
table.to_csv(pheno_out_p, index=False, sep='\t')