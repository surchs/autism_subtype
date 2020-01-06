import sys
import pathlib as pal
sys.path.append(str(pal.Path(__file__).resolve().parents[2]))
import asdfc
import argparse
import numpy as np
import patsy as pat
import pandas as pd


def abide1_nuisance(scale):
    # Hardcoded variables
    regressors = 'AGE_AT_SCAN+fd_scrubbed+SITE_ID'
    # Paths
    root_p = pal.Path(__file__).resolve().parents[2]
    sca_p = root_p / f'data/preprocessed/seed_maps/abide_2/MIST_{scale}'
    sca_t = f'sub_{{}}_ses_{{}}_rest{{}}_mist_{scale}_nocereb.npy'
    pheno_p = root_p / 'data/pheno/ABIDE2_Pheno_PSM_matched_minimum_10.tsv'
    # Output
    out_p = root_p / 'data/processed/residuals/abide_2'
    resid_p = out_p / f'abide_2_mist_{scale}_residuals_after_{"_and_".join(regressors.split("+"))}_nocereb.npy'
    if not out_p.is_dir():
        out_p.mkdir()

    pheno = pd.read_csv(pheno_p, sep='\t')
    seed_paths = [sca_p / sca_t.format(row['SUB_ID'], row['session'], row['run']) for rid, row in pheno.iterrows()]
    subject_stack = np.array([np.load(p) for p in seed_paths])
    design_matrix = pat.dmatrix(regressors, data=pheno)
    residuals = asdfc.stats.nuisance_correction(subject_stack, design_matrix, n_jobs=20)
    np.save(resid_p, residuals)
    print(f'Done saving the residuals to {resid_p}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("scale", type=int, choices=[7, 12, 20, 36, 64, 122],
                        help="specify the MIST scale the subtyping should be run on")
    args = parser.parse_args()
    abide1_nuisance(args.scale)
