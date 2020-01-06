import sys
import pathlib as pal
sys.path.append(str(pal.Path(__file__).resolve().parents[2]))
import asdfc
import argparse
import numpy as np
import pandas as pd
import patsy as pat
from tqdm import tqdm
from concurrent import futures
from sklearn import model_selection as skm


def wrap_subtype_stability(arg):
    data_stack = arg['data_stack']
    sbt_idx = arg['sbt_idx']
    dist_thr = arg['dist_thr']
    part_thr = arg['part_thr']
    regressors = arg['regressors']
    pheno = arg['pheno']
    # Regress nuisance for these individuals first
    design_matrix = pat.dmatrix(regressors, data=pheno.iloc[sbt_idx])
    residuals = asdfc.stats.nuisance_correction(data_stack[sbt_idx, ...], design_matrix, n_jobs=-1)
    # Then extract the subtype
    part, _, _ = asdfc.stats.subtype_partition(residuals, mode='core', dist_thr=dist_thr, part_thr=part_thr)
    return part


def abide1_subtype_stability_core(scale, n_cpu):
    # Hardcoded variables
    state = 1
    n_boot = 1000
    dist_thr = 0.99
    part_thr = 20
    regressors = 'AGE_AT_SCAN+fd_scrubbed+SITE_ID'
    # Paths
    root_p = pal.Path(__file__).resolve().parents[2] / 'data'
    pheno_p = root_p / 'pheno/ABIDE1_Pheno_PSM_matched_minimum_10.tsv'
    # Data
    sca_p = root_p / f'preprocessed/seed_maps/abide_1/MIST_{scale}'
    sca_t = f'sub_{{}}_ses_{{}}_run{{}}_mist_{scale}_nocereb.npy'
    # Output
    out_d = root_p / f'processed/stability/abide_1/'
    out_p = out_d / f'abide_1_subtype_stability_mist_{scale}_core_{part_thr:d}_within_{dist_thr*100:.0f}.npz'
    if not out_d.is_dir():
        out_d.mkdir()

    pheno = pd.read_csv(pheno_p, sep='\t')
    seed_paths = [sca_p / sca_t.format(row['SUB_ID'], row['session'], row['run']) for rid, row in pheno.iterrows()]
    subject_stack = np.array([np.load(p) for p in seed_paths])
    n_sub, n_vox, n_roi = subject_stack.shape
    splitter = skm.StratifiedShuffleSplit(n_splits=n_boot, test_size=0.5, random_state=state)
    asd_label = (pheno.DX_GROUP == 'Autism').values.astype(int)
    n_samples = len(asd_label)
    if not n_samples == n_sub:
        raise Exception(f'got {n_sub} subjects in residual but {n_samples} in the pheno file. Doesnt work.')
    # data_stack, mode='classic', n_subtypes=3, dist_thr=0.7, part_thr=20
    job_arg_list = [{'data_stack': subject_stack,
                     'sbt_idx': train,
                     'dist_thr': dist_thr,
                     'part_thr': part_thr,
                     'regressors': regressors,
                     'pheno': pheno}
                    for train, test in splitter.split(X=np.zeros(n_samples), y=asd_label)]
    train_indices_list, _ = zip(*list(splitter.split(X=np.zeros(n_samples), y=asd_label)))
    # decorate the subtype function
    ex = futures.ThreadPoolExecutor(max_workers=n_cpu)
    results = {run_id: res
               for run_id, res in zip(range(len(job_arg_list)),
                                      list(tqdm(ex.map(wrap_subtype_stability, job_arg_list),
                                                total=len(job_arg_list))))}

    # Store the results
    np.savez(out_p, train_idx=train_indices_list, partitions=results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("scale", type=int, choices=[7, 12, 20, 36, 64, 122],
                        help="specify the MIST scale the subtyping should be run on")
    parser.add_argument("-n_cpu", "--number_cpu", type=int, default=20,
                        help="specify the number of concurrent jobs to run")
    args = parser.parse_args()
    abide1_subtype_stability_core(args.scale, args.number_cpu)
