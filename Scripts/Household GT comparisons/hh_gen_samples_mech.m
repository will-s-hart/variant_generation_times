% Load and format samples of household generation times from the MCMC
% parameter fitting procedure.

% Need to run "Scripts/Parameter fitting/param_fit_mech.m" before this
% script.

clear all; close all; clc;

% Load output of MCMC fitting procedure

load('../../Results/param_fit_mech.mat','output')

% Generation time samples from each kept MCMC step

sample_table = output.gen_sample_table;
no_steps = size(sample_table,1);

% Infectors, infectees and generation time from each sample

infector_all = cell2mat(sample_table.infector);
infectee_all = cell2mat(sample_table.infectee);
t_gen_all = cell2mat(sample_table.t_gen);
d_gen_all = cell2mat(sample_table.d_gen);
presymp_all = cell2mat(sample_table.presymp);

% Corresponding step numbers

c1 = mat2cell((1:no_steps)',ones(no_steps,1));
step_no_all = cell2mat(cellfun(@(sample,step_no)repmat(step_no,length(sample),1),sample_table.infector,c1,'UniformOutput',false));

% Incorporate step numbers into gen_sample_all

sample_table_all = table(step_no_all,infector_all,infectee_all,t_gen_all,d_gen_all,presymp_all,'VariableNames',{'step_no','infector','infectee','t_gen','d_gen','presymp'});

% Save samples

save('../../Results/hh_gen_samples_mech.mat','sample_table','sample_table_all')