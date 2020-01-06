import sys
import pathlib as pal
sys.path.append(str(pal.Path(__file__).resolve().parents[2]))
import asdfc
import argparse
import numpy as np
import patsy as pat
import pandas as pd


def hnu1_nuisance(scale):
    # Hardcoded variables
    regressors = 'AGE_AT_SCAN+FD_scrubbed'
    # Paths
    root_p = pal.Path(__file__).resolve().parents[2]
    sca_p = root_p / f'data/preprocessed/seed_maps/hnu_1/MIST_{scale}'
    sca_t = f'sub_{{}}_run_{{}}_mist_{scale}_nocereb.npy'
    pheno_p = root_p / 'data/pheno/HNU1_full_sample_qc_passing.csv'
    # Output
    out_p = root_p / 'data/processed/residuals/hnu_1'
    resid_p = out_p / f'hnu1_session_{{}}_mist_{scale}_residuals_after_{"_and_".join(regressors.split("+"))}_nocereb.npy'
    if not out_p.is_dir():
        out_p.mkdir()

    pheno = pd.read_csv(pheno_p)
    pheno['SESSIONID'] = [int(row['SESSION'].split('_')[1]) + 1
                          for rid, row in pheno.iterrows()]

    for ses_id in range(10):
        pheno_ses = pheno.query(f'SESSIONID == {ses_id + 1}')
        seed_paths = [sca_p / sca_t.format(row['SUBID'], row['SESSIONID']) for rid, row in pheno_ses.iterrows()]
        print(f'Data for session {ses_id + 1} is : {all([s.is_file() for s in seed_paths])}')
        subject_stack = np.array([np.load(p) for p in seed_paths])
        design_matrix = pat.dmatrix(regressors, data=pheno_ses)
        residuals = asdfc.stats.nuisance_correction(subject_stack, design_matrix, n_jobs=-4)
        np.save(str(resid_p).format(ses_id + 1), residuals)
        print(f'Done saving the residuals to {str(resid_p).format(ses_id + 1)}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("scale", type=int, choices=[7, 12, 20, 36, 64, 122],
                        help="specify the MIST scale the subtyping should be run on")
    args = parser.parse_args()
    hnu1_nuisance(args.scale)