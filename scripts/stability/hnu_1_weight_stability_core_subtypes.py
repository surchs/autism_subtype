import os
import sys
import pathlib as pal
sys.path.append(str(pal.Path(__file__).resolve().parents[2]))
os.environ["MKL_NUM_THREADS"] = "2"
os.environ["NUMEXPR_NUM_THREADS"] = "2"
os.environ["OMP_NUM_THREADS"] = "2"
import asdfc
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from concurrent import futures


def hnu1_weight_stability_core(scale, n_cpu):
    # Hardcoded variables
    dist_thr = 0.99
    part_thr = 4
    regressors = 'AGE_AT_SCAN+FD_scrubbed'
    # Paths
    root_p = pal.Path(__file__).resolve().parents[2] / 'data'
    mask_root_p = root_p / 'ATLAS/MIST/'
    # Data
    resid_p = root_p / 'processed/residuals/hnu_1' / f'hnu1_session_{{}}_mist_{scale}_residuals_after_{"_and_".join(regressors.split("+"))}_nocereb.npy'
    labels_p = mask_root_p / f'Parcel_Information/MIST_{scale}_nocereb.csv'
    # Output
    out_d = root_p / 'processed/stability/hnu_1/'
    array_out_p = out_d / f'hnu_1_weight_stabilities_mist_{scale}_core_subtypes_{part_thr:d}_within_{dist_thr*100:.0f}_nocereb.npy'
    table_out_p = out_d / f'hnu_1_weight_stabilities_mist_{scale}_core_subtypes_{part_thr:d}_within_{dist_thr*100:.0f}_nocereb.tsv'
    if not out_d.is_dir():
        out_d.mkdir()

    labels = pd.read_csv(str(labels_p), sep=';')
    session_pairs = asdfc.tools.find_all_combinations(10, 2)
    residuals = np.stack([np.load(str(resid_p).format(ses+1)) for ses in range(10)], axis=-1)

    job_list = [{'data_stack': residuals,
                 'sbt_ids': sbt_ses,
                 'icc_ids': icc_ses,
                 'mode': 'core',
                 'dist_thr': dist_thr,
                 'part_thr': part_thr}
                for idx, (icc_ses, sbt_ses) in enumerate(session_pairs)]
    job_table = pd.DataFrame(data=[{'jid': jid,
                                    's_sbt': len(d['sbt_ids']),
                                    'n_icc': len(d['icc_ids']),
                                    'c_name': f'{len(d["sbt_ids"])}_{len(d["icc_ids"])}'}
                                   for jid, d in enumerate(job_list)])
    # Decorate the weight stability function
    parallel_weights = asdfc.wrappers.unpacker(asdfc.wrappers.wrap_weight_stability)
    ex = futures.ThreadPoolExecutor(max_workers=n_cpu)
    results = {run_id: res
               for run_id, res in zip(range(len(job_list)),
                                      list(tqdm(ex.map(parallel_weights, job_list), total=len(job_list))))}
    # Remap the results to array form
    res_array = np.stack([results[idx] for idx in range(len(job_list))], -1)

    print(f'Done with all {len(job_list)} jobs!')
    # We don't need to average across subtypes anymore because that happens inside the wrapper. Here we only
    # average across the sessions
    table_list = [job_table.assign(**{'icc': res_array[idx, 0, :],
                                      'wms': res_array[idx, 1, :],
                                      'bms': res_array[idx, 2, :],
                                      'n_sbt': res_array[idx, 3, :],
                                      'n_in_sbt': res_array[idx, 4, :],
                                      'n_out_sbt': res_array[idx, 5, :],
                                      'avg_size_sbt': res_array[idx, 6, :],
                                      'min_size_sbt': res_array[idx, 7, :],
                                      'network':labels.iloc[idx]['label']}) for idx in range(residuals.shape[2])]
    data_table = pd.concat(table_list)
    # Save the results
    data_table.to_csv(str(table_out_p), sep='\t')
    np.save(str(array_out_p), res_array)
    print('All Done!')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("scale", type=int, choices=[7, 12, 20, 36, 64, 122],
                        help="specify the MIST scale the subtyping should be run on")
    parser.add_argument("-n_cpu", "--number_cpu", type=int, default=20,
                        help="specify the number of concurrent jobs to run")
    args = parser.parse_args()
    hnu1_weight_stability_core(args.scale, args.number_cpu)