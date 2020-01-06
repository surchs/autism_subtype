import os
import sys
import pathlib as pal
sys.path.append(str(pal.Path(__file__).resolve().parents[2]))
os.environ["MKL_NUM_THREADS"] = "40"
os.environ["NUMEXPR_NUM_THREADS"] = "40"
os.environ["OMP_NUM_THREADS"] = "40"
import asdfc
import argparse
import pandas as pd
import nibabel as nib
from tqdm import tqdm
from concurrent import futures


def abide1_seed(scale, n_cpu):
    # Paths
    root_p = pal.Path(__file__).resolve().parents[2]
    mask_root_p = root_p / 'data/ATLAS/MIST/Parcellations/'
    mask_p = mask_root_p / 'MIST_mask_nocereb.nii.gz'
    temp_p = mask_root_p / f'MIST_{scale}_nocereb.nii.gz'
    pheno_p = root_p / 'data/pheno/ABIDE2_Pheno_QC_longitudinal_pass.tsv'
    # Data
    fc_p = root_p / 'data/preprocessed/time_series/abide_2/'
    fc_t = 'fmri_sub{:07}_session{}_rest{}_scrubbed.nii.gz'

    # Output
    out_t = f'sub_{{}}_ses_{{}}_rest{{}}_mist_{scale}_nocereb.npy'
    out_p = root_p / f'data/preprocessed/seed_maps/abide_2/MIST_{scale}'
    if not out_p.is_dir():
        out_p.mkdir()

    pheno = pd.read_csv(pheno_p, sep='\t')
    mask_i = nib.load(str(mask_p))
    atlas_i = nib.load(str(temp_p))

    path_list = [(fc_p / fc_t.format(row['SUB_ID'], row['session'], row['run']),
                  out_p / out_t.format(row['SUB_ID'], row['session'], row['run']))
                 for rid, row in pheno.iterrows()]
    for path_in, path_out in path_list:
        if not path_in.is_file():
            print(f'Path does not exist: {path_in}')
    func_in_paths, func_out_paths = zip(*path_list)

    job_list = [{'func_in_p': p_in,
                 'sca_out_p': p_out,
                 'atlas_img': atlas_i,
                 'mask_img': mask_i} for p_in, p_out in path_list]
    with futures.ThreadPoolExecutor(max_workers=n_cpu) as executor:
        futs = [executor.submit(asdfc.wrappers.wrap_seed_based_correlation, **job_args) for job_args in job_list]
        for future in tqdm(futures.as_completed(futs), total=len(job_list)):
            pass
    print('All Done!')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("scale", type=int, choices=[7, 12, 20, 36, 64, 122],
                        help="specify the MIST scale the subtyping should be run on")
    parser.add_argument("-n_cpu", "--number_cpu", type=int, default=20,
                        help="specify the number of concurrent jobs to run")
    args = parser.parse_args()
    abide1_seed(args.scale, args.number_cpu)
