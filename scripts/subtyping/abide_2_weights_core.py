import sys
import pathlib as pal
sys.path.append(str(pal.Path(__file__).resolve().parents[2]))
import asdfc
import argparse
import numpy as np
import pandas as pd


def weights_core(scale):
    # Hardcoded variables
    regressors = 'AGE_AT_SCAN+fd_scrubbed+SITE_ID'
    part_thr = 20
    dist_thr = 0.89
    # Paths
    root_p = pal.Path(__file__).resolve().parents[2] / 'data'
    label_p = root_p / 'ATLAS/MIST/Parcel_Information' / f'MIST_{scale}_nocereb.csv'
    # Data
    resid_p = root_p / 'processed/residuals/abide_2' / f'abide_2_mist_{scale}_residuals_after_{"_and_".join(regressors.split("+"))}_nocereb.npy'
    subtype_d = root_p / f'processed/subtypes/abide_1/MIST_{scale}'
    subtype_t = f'roi_{{}}_mist_{scale}_abide_1_{"_and_".join(regressors.split("+"))}_core_{part_thr:d}_within_{dist_thr*100:.0f}_nocereb'
    # Output
    out_d = root_p / f'processed/weights/abide_2/MIST_{scale}'
    if not out_d.is_dir():
        out_d.mkdir()

    label = pd.read_csv(label_p, sep=';')
    resid = np.load(resid_p)

    # Unpack the results and save them
    for scale_id in range(label.shape[0]):
        roi_number = label.iloc[scale_id]['roi']
        seed_io_name = subtype_t.format(roi_number)
        subtype_p = subtype_d / seed_io_name
        weight_out_p = out_d / f'{seed_io_name}_in_abide_2'
        subtype_arr = np.load(subtype_p.with_suffix('.npy'))
        weights = asdfc.stats.subtype_weights(resid[..., scale_id], subtype_arr)
        # Save all that
        np.save(weight_out_p.with_suffix('.npy'), weights)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("scale", type=int, choices=[7, 12, 20, 36, 64, 122],
                        help="specify the MIST scale the subtyping should be run on")
    args = parser.parse_args()
    weights_core(args.scale)
