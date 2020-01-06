#!/bin/bash
#SBATCH --account=rrg-pbellec_cpu
#SBATCH --time=1:00:00
#SBATCH --job-name=a1_pp_M20
#SBATCH --output=/project/6003287/PROJECT/abide_univariate/job_logs/%x-%j.out
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=3G

cpu=48
scale=64
# Define the set of ROIs to look at
# 122:
# MOTOR:
# roi	label
# 5	MOTnet_am
# 8	L_MOTnet_dl
# 67	MOTnet_m
# 68	MOTnet_ml
# 73	R_MOTnet_dl
# 93	MOTnet_l
# 110	MOTnet_vl # separate from the others at scale 20

# 64:
# MOTOR:
# roi	label
# 	roi	label
#   14	MOTnet_vl # Separate from the others at scale 20
#   45	MOTnet_m
#   56	MOTnet_dl



# Separate these by spaces if multiple
rois="45 56"
module load singularity/3.1
# Define some paths
data_root=/project/6003287/PROJECT/abide_univariate/
proj_root=/home/surchs/local_projects/abide_univariate
raw_root=/project/6003287/DATA/ABIDE_1/PREPROCESS_NIAK/fmri_scrubbed

# args.scale, args.confound_scale, args.confound_rois, args.number_cpu
echo Starting Partial Seed Maps
singularity run -B ${data_root},${proj_root},${raw_root} /home/surchs/singularity/abide_subtype_0_1.simg python \
    ${proj_root}/scripts/processing/abide_1_partial_seed_maps.py ${scale}  ${rois} -n_cpu ${cpu}
echo Done with Partial Seed Maps
echo Starting Partial Nuisance Regression
singularity run -B ${data_root},${proj_root},${raw_root} /home/surchs/singularity/abide_subtype_0_1.simg python \
    ${proj_root}/scripts/processing/abide_1_nuisance_regression_partial.py ${scale} ${rois}
echo Done with Partial Nuisance Regression
echo Starting Partial Subtyping
singularity run -B ${data_root},${proj_root},${raw_root} /home/surchs/singularity/abide_subtype_0_1.simg python \
    ${proj_root}/scripts/subtyping/abide_1_subtyping_core_partial.py ${scale} ${rois}
echo Done with Partial Subtyping
echo Starting Partial Weight extraction
singularity run -B ${data_root},${proj_root},${raw_root} /home/surchs/singularity/abide_subtype_0_1.simg python \
    ${proj_root}/scripts/subtyping/abide_1_weights_core_partial.py ${scale} ${rois}
echo Done with Partial Weight extraction
