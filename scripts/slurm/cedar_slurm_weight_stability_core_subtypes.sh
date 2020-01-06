#!/bin/bash
#SBATCH --account=rrg-jacquese_cpu
#SBATCH --time=1:00:00
#SBATCH --job-name=wgt_stab_core
#SBATCH --output=/project/6003287/PROJECT/abide_univariate/job_logs/%x-%j.out
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=2G

scale=20
cpu=48
# Define some paths
data_root=/project/6003287/PROJECT/abide_univariate
proj_root=/home/surchs/local_projects/abide_univariate
atlas_root=/project/6003287/ATLAS/

echo Starting Singularity
singularity run -B ${data_root},${proj_root},${atlas_root} /home/surchs/singularity/abide_subtype_0_1.simg python ${proj_root}/scripts/stability/hnu_1_weight_stability_core_subtypes.py ${scale} -n_cpu ${cpu}
echo Done with Singularity