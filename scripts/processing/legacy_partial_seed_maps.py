import os
import sys
import pathlib as pal
sys.path.append(str(pal.Path(__file__).resolve().parents[2]))
os.environ["MKL_NUM_THREADS"] = "40"
os.environ["NUMEXPR_NUM_THREADS"] = "40"
os.environ["OMP_NUM_THREADS"] = "40"
import asdfc
import argparse
import numpy as np
import pandas as pd
import nibabel as nib
from tqdm import tqdm
from concurrent import futures
from nilearn import input_data as nid


def abide1_partial_seed(scale, confound_scale, confound_rois, n_cpu):
    # Make sure the confound_rois are in an iterable
    try:
        iter(confound_rois)
    except TypeError:
        print(f'confound_rois should be an iterable. I will see if I can wrap it in a list')
        if not type(confound_rois) == int:
            raise Exception(f'confound_rois is not an iterable and also not an integer. It is {type(confound_rois)}. '
                            f'I cannot deal with this and will end here.')
        confound_rois = list(confound_rois)

    # Paths
    root_p = pal.Path(__file__).resolve().parents[2] / 'data'
    pheno_p = root_p / 'pheno/ABIDE1_Pheno_PSM_matched_minimum_10.tsv'
    atlas_p = root_p / f'ATLAS/MIST/Parcellations/MIST_{scale}_nocereb.nii.gz'
    confound_p = root_p / f'ATLAS/MIST/Parcellations/MIST_{confound_scale}_nocereb.nii.gz'
    atlas_labels_p = root_p / f'ATLAS/MIST/Parcel_Information/MIST_{scale}_nocereb.csv'
    confound_labels_p = root_p / f'ATLAS/MIST/Parcel_Information/MIST_{confound_scale}_nocereb.csv'
    mask_p = root_p / 'ATLAS/MIST/Parcellations/MIST_mask_nocereb.nii.gz'
    hier_p = root_p / 'ATLAS/MIST/Hierarchy/MIST_PARCEL_ORDER.csv'

    # Data
    fc_p = root_p / 'preprocessed/time_series/abide_1/'
    fc_t = 'fmri_sub{:07}_session{}_run{}_scrubbed.nii.gz'

    # Output
    out_t = f'sub_{{}}_ses_{{}}_run{{}}_mist_{scale}_children_of_' \
        f'{"_and_".join(map(str, confound_rois))}_mist_{confound_scale}.npy'
    out_p = root_p / f'preprocessed/seed_maps/abide_1/MIST_{scale}'
    if not out_p.is_dir():
        out_p.mkdir()

    # Load inputs
    pheno = pd.read_csv(pheno_p, sep='\t')
    mask_i = nib.load(str(mask_p))
    atlas_labels = pd.read_csv(atlas_labels_p, sep=';')
    confound_labels = pd.read_csv(confound_labels_p, sep=';')
    # Get the hierarchy
    hier = pd.read_csv(hier_p)
    ordered_labels = list()
    for i in hier[f's{scale}'].values:
        if i not in ordered_labels:
            ordered_labels.append(i)
    atlas_networks = np.array(ordered_labels)
    confound_network_list = np.array([hier.loc[hier[f's{scale}'] == i][f's{confound_scale}'].values[0]
                                      for i in ordered_labels])

    mist_mask_i = nib.load(str(mask_p))
    masker = nid.NiftiMasker(mask_img=mist_mask_i)
    masker.fit()

    # Get the atlas
    atlas_i = nib.load(str(atlas_p))
    confound_i = nib.load(str(confound_p))
    atlas = masker.transform(atlas_i)
    confound = masker.transform(confound_i)

    # Find the regions and give their names
    # Pick the atlas child networks that belong to the confound regions
    confound_names = list(confound_labels.query('roi in @confound_rois')['label'].values)
    children = np.array([atlas_roi
                         for confound_roi in confound_rois
                         for atlas_roi in atlas_networks[confound_network_list == confound_roi]])
    children_names = list(atlas_labels.query('roi in @children')['label'].values)
    print(f'I am generating seed maps for the children of {confound_names} from MIST_{confound_scale} '
          f'in MIST_{scale}. These children are:\n{children_names}')
    # Build the new atlases
    confound_mask = np.array([1 if confound[0, idx] in confound_rois else 0 for idx in range(confound.shape[1])])
    subset_atlas = np.array([0 if not atlas[0, idx] in children else atlas[0, idx] for idx in range(atlas.shape[1])])
    confound_mask_i = masker.inverse_transform(confound_mask)
    subset_atlas_i = masker.inverse_transform(subset_atlas)

    path_list = [(fc_p / fc_t.format(row['SUB_ID'], row['session'], row['run']),
                  out_p / out_t.format(row['SUB_ID'], row['session'], row['run']))
                 for rid, row in pheno.iterrows()]
    for path_in, path_out in path_list:
        if not path_in.is_file():
            print(f'Path does not exist: {path_in}')
    func_in_paths, func_out_paths = zip(*path_list)

    job_list = [{'func_in_p': p_in,
                 'sca_out_p': p_out,
                 'atlas_img': subset_atlas_i,
                 'mask_img': mask_i,
                 'confound_img': confound_mask_i} for p_in, p_out in path_list]
    with futures.ThreadPoolExecutor(max_workers=n_cpu) as executor:
        futs = [executor.submit(asdfc.wrappers.wrap_seed_based_correlation, **job_args) for job_args in job_list]
        for future in tqdm(futures.as_completed(futs), total=len(job_list)):
            pass
    print('All Done!')


if __name__ == "__main__":
    # confound_scale, confound_rois
    parser = argparse.ArgumentParser()
    parser.add_argument("scale", type=int, choices=[7, 12, 20, 36, 64, 122],
                        help="specify the MIST scale the subtyping should be run on")
    parser.add_argument("confound_scale", type=int, choices=[7, 12, 20, 36, 64, 122],
                        help="specify the MIST scale the subtyping regressors should come from")
    parser.add_argument('-croi', '--confound_rois', nargs='+', type=int)
    parser.add_argument("-n_cpu", "--number_cpu", type=int, default=20,
                        help="specify the number of concurrent jobs to run")
    args = parser.parse_args()
    abide1_partial_seed(args.scale, args.confound_scale, args.confound_rois, args.number_cpu)

