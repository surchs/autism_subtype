%%% NKI TRT preprocessing pipeline
% Script to run a preprocessing pipeline analysis on the HCP database.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting input/output files %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
root_p = '/project/6003287/DATA/ABIDE_1';
site_name = 'site';
load([root_p filesep 'Support_NIAK' filesep site_name '_niak.mat']);
opt.flag_test = false;
% [pipeline,opt] = niak_pipeline_fmri_preprocess(files_in,opt);
