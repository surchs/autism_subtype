import sys
import pathlib as pal
sys.path.append(str(pal.Path(__file__).resolve().parents[2]))
import json
import asdfc
import argparse
import numpy as np
import nibabel as nib
from nilearn import input_data


def subtyping_fixed(scale, n_sbt):
    # Hardcoded variables
    regressors = 'AGE_AT_SCAN+fd_scrubbed+SITE_ID'
    # Paths
    root_p = pal.Path(__file__).resolve().parents[2] / 'data'
    mask_p = root_p / 'ATLAS/MIST/Parcellations/MIST_mask_nocereb.nii.gz'
    temp_p = root_p / f'ATLAS/MIST/Parcellations/MIST_{scale}_nocereb.nii.gz'
    # Data
    resid_p = root_p / 'processed/residuals/abide_1' / f'abide_1_mist_{scale}_residuals_after_{"_and_".join(regressors.split("+"))}_nocereb.npy'
    # Output
    out_d = root_p / f'processed/subtypes/abide_1/MIST_{scale}'
    out_t = f'roi_{{}}_mist_{scale}_abide_1_{"_and_".join(regressors.split("+"))}_fixed_{n_sbt}_nocereb'
    if not out_d.is_dir():
        out_d.mkdir()

    resid = np.load(resid_p)
    mask_i = nib.load(str(mask_p))
    voxel_masker = input_data.NiftiMasker(mask_img=mask_i, standardize=False)
    voxel_masker.fit()

    # Compute the subtypes
    partition, distance, order = asdfc.stats.subtype_partition(resid, mode='classic', n_subtypes=n_sbt)
    subtype_maps = asdfc.stats.subtype_maps(resid, partition)

    # Unpack the results and save them
    for scale_id in range(scale):
        scale_out_p = out_d / out_t.format(scale_id+1)
        json_d = {'part': list(map(int, partition[:, scale_id])),
                  'order': list(map(int, order[:, scale_id])),
                  'subtype': 'fixed',
                  'n_subtypes': n_sbt,
                  'residual_used': str(resid_p),
                  'atlas_used': str(temp_p)}
        subtype_arr = subtype_maps[scale_id]
        subtypes_i = voxel_masker.inverse_transform(subtype_arr)
        # Save all that
        with scale_out_p.with_suffix('.json').open('w') as f:
            json.dump(json_d, f)
        np.save(scale_out_p.with_suffix('.npy'), subtype_arr)
        nib.save(subtypes_i, str(scale_out_p.with_suffix('.nii.gz')))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("scale", type=int, choices=[7, 12, 20, 36, 64, 122],
                        help="specify the MIST scale the subtyping should be run on")
    parser.add_argument("-n_sbt", "--n_subtype", type=int, default=5,
                        help="specify the number of subtypes that should be extracted for each network")
    args = parser.parse_args()
    subtyping_fixed(args.scale, args.n_subtype)
