#!/bin/bash
#SBATCH --account=rrg-pbellec_cpu
#SBATCH --time=1:00:00
#SBATCH --job-name=abide2_scrub
#SBATCH --output=/project/6003287/PROJECT/abide_univariate/job_logs/%x-%j.out
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=2G

cpu=30
module load singularity/3.2
# Define some paths
data_root=/project/6003287/DATA/ABIDE_2/PREPROCESS_NIAK
proj_root=/home/surchs/local_projects/abide_univariate

echo Starting Scrubbing
singularity run -B ${data_root},${proj_root} /home/surchs/singularity/abide_subtype_0_1.simg python ${proj_root}/scripts/preprocessing/abide_2_niak_scrubbing.py -n_cpu ${cpu}
echo Done with Scrubbing
