#!/bin/bash
#SBATCH --account=rpp-aevans-ab_cpu
#SBATCH --time=1:00:00
#SBATCH --job-name=a1_proc_M122
#SBATCH --output=/lustre03/project/6008063/surchs/PROJECTS/abide_univariate/job_logs/%x-%j.out
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=3G

scale=20
cpu=48
module load singularity/3.1
# Define some paths
data_root=/lustre03/project/6008063/surchs/PROJECTS/abide_univariate
proj_root=/home/surchs/local_projects/abide_univariate
raw_root=/lustre03/project/6008063/surchs/DATA/ABIDE_1/PREPROCESS_NIAK/fmri_scrubbed

echo Starting Seed Maps
#singularity run -B ${data_root},${proj_root},${raw_root} /home/surchs/singularity/abide_subtype_0_1.simg python ${proj_root}/scripts/processing/abide_1_seed_maps.py ${scale} -n_cpu ${cpu}
echo Done with Seed Maps
echo Starting Nuisance Regression
singularity run -B ${data_root},${proj_root},${raw_root} /home/surchs/singularity/abide_subtype_0_1.simg python ${proj_root}/scripts/processing/abide_1_nuisance_regression.py ${scale}
echo Done with Nuisance Regression
echo Starting Subtyping
singularity run -B ${data_root},${proj_root},${raw_root} /home/surchs/singularity/abide_subtype_0_1.simg python ${proj_root}/scripts/subtyping/abide_1_subtyping_core.py ${scale}
echo Done with Subtyping
echo Starting Weight extraction
singularity run -B ${data_root},${proj_root},${raw_root} /home/surchs/singularity/abide_subtype_0_1.simg python ${proj_root}/scripts/subtyping/abide_1_weights_core.py ${scale}
echo Done with Weight extraction
