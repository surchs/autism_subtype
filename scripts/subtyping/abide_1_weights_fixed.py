import sys
import pathlib as pal
sys.path.append(str(pal.Path(__file__).resolve().parents[2]))
import asdfc
import argparse
import numpy as np
import nibabel as nib
from nilearn import input_data


def weights_fixed(scale, n_sbt):
    # Hardcoded variables
    regressors = 'AGE_AT_SCAN+fd_scrubbed+SITE_ID'
    # Paths
    root_p = pal.Path(__file__).resolve().parents[2] / 'data'
    mask_p = root_p / 'ATLAS/MIST/Parcellations/MIST_mask_nocereb.nii.gz'
    # Data
    resid_p = root_p / 'processed/residuals/abide_1' / f'abide_1_mist_{scale}_residuals_after_{"_and_".join(regressors.split("+"))}_nocereb.npy'
    subtype_d = root_p / f'processed/subtypes/abide_1/MIST_{scale}'
    subtype_t = f'roi_{{}}_mist_{scale}_abide_1_{"_and_".join(regressors.split("+"))}_fixed_{n_sbt}_nocereb'
    # Output
    out_d = root_p / f'processed/weights/abide_1/MIST_{scale}'
    if not out_d.is_dir():
        out_d.mkdir()

    resid = np.load(resid_p)
    mask_i = nib.load(str(mask_p))
    voxel_masker = input_data.NiftiMasker(mask_img=mask_i, standardize=False)
    voxel_masker.fit()

    # Compute the subtypes
    partition, distance, order = asdfc.stats.subtype_partition(resid, mode='classic', n_subtypes=5)
    subtype_maps = asdfc.stats.subtype_maps(resid, partition)

    # Unpack the results and save them
    for scale_id in range(len(subtype_maps)):
        seed_io_name = subtype_t.format(scale_id+1)
        subtype_p = subtype_d / seed_io_name
        weight_out_p = out_d / seed_io_name
        subtype_arr = np.load(subtype_p.with_suffix('.npy'))
        weights = asdfc.stats.subtype_weights(resid[..., scale_id], subtype_arr)
        # Save all that
        np.save(weight_out_p.with_suffix('.npy'), weights)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("scale", type=int, choices=[7, 12, 20, 36, 64, 122],
                        help="specify the MIST scale the subtyping should be run on")
    parser.add_argument("-n_subtype", "--n_sbt", type=int, default=5,
                        help="specify the number of subtypes that should be extracted for each network")
    args = parser.parse_args()
    weights_fixed(args.scale)

