
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting input/output files %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
root_p = '/home/surchs/Projects/abide_univariate/data/pheno/';
load([root_p filesep 'abide1_case_control_support.mat']);
opt.flag_test = false;
[pipeline,opt] = niak_pipeline_glm_connectome(files_in,opt);
