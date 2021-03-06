import os
import sys
import pathlib as pal
sys.path.append(str(pal.Path(__file__).resolve().parents[2]))
os.environ["MKL_NUM_THREADS"] = "40"
os.environ["NUMEXPR_NUM_THREADS"] = "40"
os.environ["OMP_NUM_THREADS"] = "40"
import asdfc
import argparse
from tqdm import tqdm
from concurrent import futures


def cnv_scrub(n_cpu):
    # Paths
    project_root_p = pal.Path('/mnt/data_sq/cisl/paper_16p22q/preprocessed_data/time_series/')
    prep_d = project_root_p / '16p_time_series/'
    out_d = project_root_p / '16p_time_series/fmri_scrubbed/'
    if not out_d.is_dir():
        out_d.mkdir()

    # Find all the images in the sub-sessions
    img_paths = list(prep_d.glob('**/fmri/*.nii.gz'))
    img_names = [p.name.split('.')[0] for p in img_paths]
    extra_paths, extra_available = zip(
        *[(img_path.parent / f'{img_name}_extra.mat', (img_path.parent / f'{img_name}_extra.mat').exists())
          for img_name, img_path in zip(img_names, img_paths)])
    out_paths = [out_d / f'{img_name}_scrubbed.nii.gz' for img_name in img_names]
    print(f'{"" if all(extra_available) else "not"} all extra.mat files are available for {len(img_names)} files found.')

    # Make a list of job dictionaries
    job_list = [{'img_p': i, 'extra_p': e, 'out_p': o, 'clobber': False} for i, e, o in
                zip(img_paths, extra_paths, out_paths)]
    with futures.ThreadPoolExecutor(max_workers=n_cpu) as executor:
        futs = [executor.submit(asdfc.data.niak_scrubbing, **job_args) for job_args in job_list]
        for future in tqdm(futures.as_completed(futs), total=len(job_list)):
            pass
    print(f'Finished the scrubbing: all jobs succeeded.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n_cpu", "--number_cpu", type=int, default=20,
                        help="specify the number of concurrent jobs to run")
    args = parser.parse_args()
    cnv_scrub(args.number_cpu)
