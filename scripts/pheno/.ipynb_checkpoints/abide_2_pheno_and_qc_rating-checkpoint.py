import re
import json
import pandas as pd
import pathlib as pal

root_p = pal.Path('/home/surchs/temp/abide_univariate_data')
pheno_p = root_p / 'pheno/ABIDE2_Pheno_prepared_motion.tsv'
pheno_qc_p = root_p / 'pheno/ABIDE2_Pheno_QC.tsv'
pheno_qc_pass_p = root_p / 'pheno/ABIDE2_Pheno_QC_pass.tsv'
report_d = pal.Path('/home/surchs/ABIDE_2_QC_Reports')

pheno = pd.read_csv(pheno_p, sep='\t')
reports = list(report_d.glob('*.json'))

pheno['QC'] = 'Unrated'
pheno['QC_comment'] = ''
for report_p in reports:
    with report_p.open('r') as f:
        qc = json.load(f)
    for sub in qc['subjects']:
        sub_id = re.search(r'(?<=sub00)\d+', sub).group()
        sub_idx = pheno.query(f'SUB_ID=={sub_id}').index
        rating = qc['subjects'][sub]['status']
        comment = qc['subjects'][sub]['comments']
        if rating == '':
            print('Something Something no Rating')
            print(sub_id, pheno.query(f'SUB_ID=={sub_id}')['SITE_ID'].values[0], report_p.name)
        pheno.loc[sub_idx, ['QC', 'QC_comment']] = rating, comment

pheno.to_csv(pheno_qc_p, index=False, sep='\t')
# Make conditions for the subjects that are considered to pass QC
available_pheno = pheno.query('QC=="Pass" '
                              'and frames_ok > 40 '
                              'and fd_scrubbed < 0.5 '
                              'and SEX == "male" '
                              'and AGE_AT_SCAN < 40')
# Remove duplicates by picking the run with the lowest motion
dup_sub = available_pheno[available_pheno.duplicated(['SUB_ID'])].SUB_ID.unique()
for sub in dup_sub:
    a = available_pheno .query(f'SUB_ID=={sub}')
    min_idx = a.fd_scrubbed.idxmin()
    drop_idx = [i for i in a.index if not i == min_idx]
    available_pheno .drop(index=drop_idx, inplace=True)

# Check that at least 5 of both categories exist in each site after the other qc criteria
good_sites = [s for s in available_pheno.SITE_ID.unique()
              if available_pheno.query(f'SITE_ID=="{s}" and DX_GROUP=="Autism"').shape[0] > 5
              and available_pheno.query(f'SITE_ID=="{s}" and DX_GROUP=="Control"').shape[0] > 5]
available_site_pheno = available_pheno.query('SITE_ID in @good_sites')

available_site_pheno.reset_index(inplace=True)
available_site_pheno.to_csv(pheno_qc_pass_p, index=False, sep='\t')
