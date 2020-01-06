import os
import sys
import pathlib as pal
sys.path.append(str(pal.Path(__file__).resolve().parents[2]))
os.environ["MKL_NUM_THREADS"] = "40"
os.environ["NUMEXPR_NUM_THREADS"] = "40"
os.environ["OMP_NUM_THREADS"] = "40"
import re
import asdfc
import argparse
import nibabel as nib
from tqdm import tqdm
from concurrent import futures


def cnv_seed(scale, n_cpu):
    # Paths
    root_p = pal.Path(__file__).resolve().parents[2]
    mask_root_p = root_p / 'data/ATLAS/MIST/Parcellations/'
    mask_p = mask_root_p / 'MIST_mask_nocereb.nii.gz'
    temp_p = mask_root_p / f'MIST_{scale}_nocereb.nii.gz'
    fc_p = pal.Path(
        '/mnt/data_sq/cisl/paper_16p22q/preprocessed_data/time_series/16p_time_series/fmri_scrubbed')
    out_p = root_p / f'data/preprocessed/seed_maps/cnv16p/MIST_{scale}/'

    # Output
    if not out_p.is_dir():
        out_p.mkdir()

    mask_i = nib.load(str(mask_p))
    atlas_i = nib.load(str(temp_p))

    func_in_paths = list(fc_p.glob('*.nii.gz'))
    func_out_paths = [out_p / '{}_mist_{}_nocereb.npy'.format(re.search(r'(?<=fmri_)s\S+(?=_rest_scrubbed.nii.gz)',
                                                                        f.name).group(), scale) for f in func_in_paths]

    for path_in, path_out in zip(func_in_paths, func_out_paths):
        if not path_in.is_file():
            print(f'Path does not exist: {path_in}')

    job_list = [{'func_in_p': p_in,
                 'sca_out_p': p_out,
                 'atlas_img': atlas_i,
                 'mask_img': mask_i} for p_in, p_out in zip(func_in_paths, func_out_paths)]

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
    cnv_seed(args.scale, args.number_cpu)
