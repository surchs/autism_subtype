import pathlib as pal
import sys
import pathlib as pal
sys.path.append(str(pal.Path(__file__).resolve().parents[2]))
import json
import asdfc
import argparse
import numpy as np
import pandas as pd
import nibabel as nib
from nilearn import input_data


def subtyping_core(scale=20):
    # Hardcoded variables
    regressors = 'AGE_AT_SCAN+fd_scrubbed+SITE_ID'
    part_thr = 40
    dist_thr = 1

    # Paths
    root_p = pal.Path(__file__).resolve().parents[2] / 'data'
    mask_p = root_p / 'ATLAS/MIST/Parcellations/MIST_mask_nocereb.nii.gz'
    temp_p = root_p / f'ATLAS/MIST/Parcellations/MIST_{scale}_nocereb.nii.gz'
    label_p = root_p / 'ATLAS/MIST/Parcel_Information' / f'MIST_{scale}_nocereb.csv'
    # Data
    resid_p = root_p / 'processed/residuals/abide_2' / f'abide_2_mist_{scale}_residuals_after_{"_and_".join(regressors.split("+"))}_nocereb.npy'
    # Output
    out_d = root_p / f'processed/subtypes/abide_2/MIST_{scale}'
    out_t = f'roi_{{}}_mist_{scale}_abide_2_{"_and_".join(regressors.split("+"))}_core_{part_thr:d}_within_{dist_thr*100:.0f}_nocereb'
    if not out_d.is_dir():
        out_d.mkdir()

    label = pd.read_csv(label_p, sep=';')
    resid = np.load(resid_p)
    mask_i = nib.load(str(mask_p))
    voxel_masker = input_data.NiftiMasker(mask_img=mask_i, standardize=False)
    voxel_masker.fit()

    # Compute the subtypes
    print('Beginning to Generate the Subtypes')
    partition, distance, order = asdfc.stats.subtype_partition(resid, mode='core', dist_thr=dist_thr, part_thr=part_thr)
    subtype_maps = asdfc.stats.subtype_maps(resid, partition)

    # Unpack the results and save them
    for scale_id in range(label.shape[0]):
        roi_number = label.iloc[scale_id]['roi']
        scale_out_p = out_d / out_t.format(roi_number)
        json_d = {'part': list(map(int, partition[:, scale_id])),
                  'order': list(map(int, order[:, scale_id])),
                  'subtype': 'core',
                  'dist_thr': float(dist_thr),
                  'part_thr': int(part_thr),
                  'residual_used': str(resid_p),
                  'atlas_used': str(temp_p)}
        subtype_arr = subtype_maps[scale_id]
        subtypes_i = voxel_masker.inverse_transform(subtype_arr)
        # Save all that
        with scale_out_p.with_suffix('.json').open('w') as f:
            json.dump(json_d, f)
        np.save(scale_out_p.with_suffix('.npy'), subtype_arr)
        nib.save(subtypes_i, str(scale_out_p.with_suffix('.nii.gz')))
    print('Done with the subtyping')


if __name__ == "__main__":
    print(pal.Path(__file__).resolve().parents[2])
    parser = argparse.ArgumentParser()
    parser.add_argument("scale", type=int, choices=[7, 12, 20, 36, 64, 122],
                        help="specify the MIST scale the subtyping should be run on")
    args = parser.parse_args()
    subtyping_core(args.scale)
