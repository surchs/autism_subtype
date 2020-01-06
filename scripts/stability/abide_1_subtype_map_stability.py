import sys
import pathlib as pal
sys.path.append(str(pal.Path(__file__).resolve().parents[2]))
import tqdm
import asdfc
import numpy as np
import pandas as pd
import pathlib as pal
import nibabel as nib
import patsy as pat
from matplotlib import gridspec
from scipy import cluster as scl
from matplotlib import pyplot as plt
from nilearn import input_data as nii
from nilearn import plotting as nlp
import random
import itertools as it
from concurrent import futures
from sklearn import metrics as skm


# Paths
scale = 20
regressors = 'AGE_AT_SCAN+fd_scrubbed+SITE_ID'
n_fold_combos = 1000
dist_thr = 0.99
part_thr = 20

root_p = pal.Path(__file__).resolve().parents[2] / 'data'
pheno_p = root_p / 'pheno/ABIDE1_Pheno_PSM_matched_minimum_10.tsv'
mask_p = root_p / 'ATLAS/MIST/Parcellations/MIST_mask_nocereb.nii.gz'
stab_p = root_p / 'processed/stability/abide_I don't know how to 1/abide_1_subtype_stability_mist_20_core_20_within_99.npz'
sca_p = root_p / f'preprocessed/seed_maps/abide_1/MIST_{scale}'
sca_t = f'sub_{{}}_ses_{{}}_run{{}}_mist_{scale}_nocereb.npy'
# Output
out_d = root_p / f'processed/stability/abide_1/'
out_p = out_d / f'abide_1_subtype_map_stability_mist_{scale}_core_{part_thr:d}_within_{dist_thr*100:.0f}.npy'
if not out_d.is_dir():
    out_d.mkdir()

mask_i = nib.load(str(mask_p))
masker = nii.NiftiMasker(mask_img=mask_i)
masker.fit()
stab = np.load(stab_p)
pheno = pd.read_csv(pheno_p, sep='\t')
seed_paths = [sca_p / sca_t.format(row['SUB_ID'], row['session'], row['run']) for rid, row in pheno.iterrows()]
subject_stack = np.array([np.load(p) for p in seed_paths])

partitions = stab['partitions'][()]
train_indices = stab['train_idx']

n_samples = pheno.shape[0]
n_boot = len(partitions.keys())
n_networks = subject_stack.shape[2]

max_sim_mat = np.zeros(shape=(n_fold_combos, n_networks, 2))
fold_combinations = list(it.combinations(range(n_boot), 2))
# A random subset of these combinations
rand_fold_comb = random.sample(fold_combinations,n_fold_combos)
total_jobs = n_networks * n_fold_combos
with tqdm.tqdm(total=total_jobs) as pbar:
    for nid in range(n_networks):
        for fid, (f1, f2) in enumerate(rand_fold_comb):
            t1 = train_indices[f1]
            t2 = train_indices[f2]
            p1 = partitions[f1][:, nid]
            p2 = partitions[f2][:, nid]
            dm1 = pat.dmatrix(regressors, data=pheno.iloc[t1])
            r1 = asdfc.stats.nuisance_correction(subject_stack[t1, :, nid], dm1, n_jobs=-1)
            dm2 = pat.dmatrix(regressors, data=pheno.iloc[t2])
            r2 = asdfc.stats.nuisance_correction(subject_stack[t2, :, nid], dm2, n_jobs=-1)

            # Compute the seed maps
            m1 = np.stack([np.mean(r1[p1==i], 0) for i in np.unique(p1[p1!=0])], -1)
            m2 = np.stack([np.mean(r2[p2==i], 0) for i in np.unique(p2[p2!=0])], -1)
            # Correlate the two maps
            c = asdfc.stats.corr2_coeff(m1.T, m2.T)
            # Find maximal correlation of f1 subtypes with f2 subtypes
            max_c = np.max(c, 1)
            # Average them
            max_sim_mat[fid, nid, 0] = np.mean(max_c)
            # And get the SD
            max_sim_mat[fid, nid, 1] = np.std(max_c)
            pbar.update(1)

# Save the output
np.save(out_p, max_sim_mat)