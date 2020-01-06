#!/bin/bash
#SBATCH --account=rpp-aevans-ab_cpu
#SBATCH --time=1:00:00
#SBATCH --job-name=sbt_stab_fix
#SBATCH --output=/project/6003287/PROJECT/abide_univariate/job_logs/%x-%j.out
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=1G

scale=12
cpu=48
module load singularity/3.1
# Define some paths
data_root=/project/6003287/PROJECT/abide_univariate
proj_root=/home/surchs/local_projects/abide_univariate

echo Starting Singularity
singularity run -B ${data_root},${proj_root} /home/surchs/singularity/abide_subtype_0_1.simg python ${proj_root}/scripts/stability/abide_1_subtype_stability_fixed.py ${scale} -n_cpu ${cpu}
echo Done with Singularity