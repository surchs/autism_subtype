import sys
import pathlib as pal
sys.path.append(str(pal.Path(__file__).resolve().parents[2]))
import asdfc
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from concurrent import futures
from sklearn import model_selection as skm


def wrap_subtype_stability_fixed(arg):
    data_stack = arg['data_stack']
    sbt_idx = arg['sbt_idx']
    n_sbt = arg['n_subtypes']
    part, _, _ = asdfc.stats.subtype_partition(data_stack[sbt_idx, ...], mode='classic', n_subtypes=n_sbt)
    return part


def abide1_subtype_stability_fixed(scale, n_cpu):
    # Hardcoded variables
    state = 1
    n_boot = 1000
    n_subtypes = 3
    regressors = 'AGE_AT_SCAN+FD_scrubbed+SITE_ID'

    # Paths
    root_p = pal.Path(__file__).resolve().parents[2] / 'data'
    pheno_p = root_p / 'pheno/abide_consensus_model_ados.csv'
    # Data
    resid_p = root_p / f'processed/residuals/abide_1/abide_1_mist_{scale}_residuals_after_{"_and_".join(regressors.split("+"))}.npy'
    # Output
    out_d = root_p / f'processed/stability/abide_1/'
    out_p = out_d / f'abide_1_subtype_stability_mist_{scale}_fixed_{n_subtypes}_subtype.npy'
    info_out_p = out_d / f'abide_1_subtype_stability_mist_{scale}_fixed_{n_subtypes}_subtype_info.npy'
    if not out_d.is_dir():
        out_d.mkdir()

    pheno = pd.read_csv(pheno_p)
    resid = np.load(resid_p)
    splitter = skm.StratifiedShuffleSplit(n_splits=n_boot, test_size=0.5, random_state=state)
    asd_label = (pheno.DX_GROUP == 'Autism').values.astype(int)
    n_samples = len(asd_label)
    # data_stack, mode='classic', n_subtypes=3, dist_thr=0.7, part_thr=20
    job_arg_list = [{'data_stack': resid,
                     'sbt_idx': train,
                     'n_subtypes': 3}
                    for train, test in splitter.split(X=np.zeros(n_samples), y=asd_label)]
    train_indices_list, _ = zip(*list(splitter.split(X=np.zeros(n_samples), y=asd_label)))
    # decorate the subtype function
    ex = futures.ThreadPoolExecutor(max_workers=n_cpu)
    results = {run_id: res
               for run_id, res in zip(range(len(job_arg_list)),
                                      list(tqdm(ex.map(wrap_subtype_stability_fixed, job_arg_list), total=len(job_arg_list))))}
    adjacency_mat = np.zeros(shape=(n_samples, n_samples, scale))
    count_mat = np.zeros(shape=(n_samples, n_samples, scale))
    # Dimensions: n_boot, seeds, len([n_sub_unassigned, n_sub_assigned, n_subtype])
    info_mat = np.zeros(shape=(n_boot, scale, 3))
    for rid in range(len(job_arg_list)):
        partition = results[rid]

        adjacency_stack = np.array([[partition[:, i] == p
                                     for p in partition[:, i]]
                                    for i in range(scale)])

        n_empty = np.array([np.sum(partition[:, scale_id] == 0) for scale_id in range(partition.shape[1])])
        n_assign = np.array([np.sum(partition[:, scale_id] != 0) for scale_id in range(partition.shape[1])])
        n_subtypes = np.array([len(np.unique(partition[:, scale_id])) for scale_id in range(partition.shape[1])])
        info_mat[rid, ...] = np.stack([n_empty, n_assign, n_subtypes], -1)

        train_ind = train_indices_list[rid]
        col_train, row_train = np.ix_(train_ind, train_ind)
        for scale_id in range(scale):
            count_mat[col_train, row_train, scale_id] += 1
            adjacency_mat[col_train, row_train, scale_id] += adjacency_stack[scale_id, ...]

    average_mat = adjacency_mat / count_mat
    average_mat[np.isnan(average_mat)] = 0

    np.save(str(out_p), average_mat)
    np.save(str(info_out_p), info_mat)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("scale", type=int, choices=[7, 12, 20, 36, 64, 122],
                        help="specify the MIST scale the subtyping should be run on")
    parser.add_argument("-n_cpu", "--number_cpu", type=int, default=20,
                        help="specify the number of concurrent jobs to run")
    args = parser.parse_args()
    abide1_subtype_stability_fixed(args.scale, args.number_cpu)