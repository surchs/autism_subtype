import re
import pandas as pd
import pathlib as pal

root_p = pal.Path('/home/surchs/data/surchs/ABIDE/ABIDE_1/')
raw_pheno_p = root_p / 'Pheno/Phenotypic_V1_0b.csv'
prepared_pheno_p = '../../data/pheno/ABIDE1_Pheno_prepared.tsv'

raw_pheno = pd.read_csv(raw_pheno_p)

"""
Fix some disagreements between the data on disk and the pheno folder
- CMU_a and CMU_b aren't separated in the pheno file but have different acquisition parameters
- MaxMun_a - MaxMun_d aren't separated in the pheno file but have different acquisition parameters
- CALTECH, OLIN, STANFORD, TRINITY, YALE, MAXMUN, LEUVEN, PITT are camelcase on disk but allcaps in the pheno file
- There are guys in the FILE_ID column that have `no_filename` there. Need to fix this so I can reliably use that column
  later on

I fix this by:
- adding a new column in the pheno file called `site_name` that is camelcase for the sites that need it
- fix the missing or incorrect entries in the FILE_ID column by checking the site name in the raw data
- for subjects in `MAX_MUN` and `CMU`, the subscript assignment is encoded in the `FILE_ID` column. So I will just parse 
  that and assign the `site_name` column based on that
"""
# Replace the numerical values by the explicit string values
raw_pheno.DX_GROUP.replace({1: 'Autism', 2: 'Control'}, inplace=True)
raw_pheno.DSM_IV_TR.replace({0: 'Control', 1: 'Autism', 2: 'Aspergers', 3:'PDD-NOS', 4:'Aspergers or PDD-NOS'}, inplace=True)
raw_pheno.SEX.replace({1: 'Male', 2: 'Female'}, inplace=True)
raw_pheno.ADI_R_RSRCH_RELIABLE.replace({0: 'not research reliable', 1: 'research reliable'}, inplace=True)
raw_pheno.ADOS_RSRCH_RELIABLE.replace({0: 'not research reliable', 1: 'research reliable'}, inplace=True)
raw_pheno.EYE_STATUS_AT_SCAN.replace({1: 'open', 2: 'closed'}, inplace=True)

# Replace missing values with NaN
raw_pheno.replace({-9999: pd.np.nan}, inplace=True)

# Camelcase sites that need it
raw_pheno['site_name'] = raw_pheno.SITE_ID.values
rename_sites = ['CALTECH', 'LEUVEN_1', 'LEUVEN_2', 'OLIN', 'PITT', 'STANFORD', 'TRINITY', 'YALE']
raw_pheno.site_name.replace({old_name:old_name.title() for old_name in rename_sites}, inplace=True)

# Fix the FILE_ID column
search_p = root_p / 'RAW'
for rid, row in raw_pheno.iterrows():
    file_id = re.match(r'([A-Za-z]+_*[a-z]*(?=_[0-9]*))', row.FILE_ID).group()
    if 'no' in file_id:
        glob_results = list(search_p.glob('*/sub-{:07}'.format(row.SUB_ID)))
        if glob_results:
            site = glob_results[0].parent.name
            file_id = '{}_{:07}'.format(site, row.SUB_ID)
            # Write it back into the pheno
            raw_pheno.loc[rid, 'FILE_ID'] = file_id
        else:
            print('Subject {} from {} is messed up'.format(row.SUB_ID, row.SITE_ID))

# subscript sites that need it with _a and _b ...
fix_index = (raw_pheno.SITE_ID == 'CMU') | (raw_pheno.SITE_ID == 'MAX_MUN')
raw_pheno.loc[fix_index, 'site_name'] = raw_pheno.loc[fix_index, 'FILE_ID'].str.extract(r'([A-Za-z]+_[a-z]*)').values

# Save the pheno file
raw_pheno.to_csv(str(prepared_pheno_p), index=False, sep='\t')