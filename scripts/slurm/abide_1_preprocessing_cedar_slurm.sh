#!/bin/bash
#SBATCH --account=rpp-aevans-ab_cpu
#SBATCH --time=0:30:00
#SBATCH --job-name=cnv_scrub
#SBATCH --output=/project/6003287/PROJECT/abide_univariate/job_logs/%x-%j.out
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=1G

cpu=48
module load singularity/3.1
# Define some paths
data_root=/project/6003287/DATA/ABIDE_1/PREPROCESS_NIAK
proj_root=/home/surchs/local_projects/abide_univariate

echo Starting Scrubbing
singularity run -B ${data_root},${proj_root} /home/surchs/singularity/abide_subtype_0_1.simg python ${proj_root}/scripts/preprocessing/abide_1_niak_scrubbing.py -n_cpu ${cpu}
echo Done with Scrubbing
