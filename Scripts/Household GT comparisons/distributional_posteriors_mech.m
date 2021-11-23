% Obtain discrete- and continuous-time posterior distributions of the
% household generation time (rather than just the mean and standard
% deviation of that distribution as considered elsewhere) for the Alpha and
% Delta variants.

% Need to run "Scripts/Parameter fitting/hh_gen_samples_mech.m" before this
% script.

clear all; close all; clc;

% Load data

load('../../Data/data.mat','data_struct_observed')

% Load output of MCMC fitting procedure

load('../../Results/hh_gen_samples_mech.mat','sample_table','sample_table_all')

% Vector of edges used to determine number of GT samples in different
% categories per step of the MCMC chain

no_steps = size(sample_table,1);
edges = (0.5:(no_steps+1))';

% GT samples by variant

variant = data_struct_observed.variant;

Alpha = (variant==1);
Delta = (variant==2);

Alpha_inds = Alpha(sample_table_all.infector);
Delta_inds = Delta(sample_table_all.infector);

sample_table_Alpha = sample_table_all(Alpha_inds,:);
sample_table_Delta = sample_table_all(Delta_inds,:);

[step_freqs_Alpha,~] = histcounts(sample_table_Alpha.step_no,edges);
step_freqs_Alpha = step_freqs_Alpha';

times_all_Alpha = sample_table_Alpha.t_gen; %all continuous-time samples
days_cell_Alpha = mat2cell(sample_table_Alpha.d_gen,step_freqs_Alpha); %cell array of discrete-time samples from each MCMC step

[step_freqs_Delta,~] = histcounts(sample_table_Delta.step_no,edges);
step_freqs_Delta = step_freqs_Delta';

times_all_Delta = sample_table_Delta.t_gen; %all continuous-time samples
days_cell_Delta = mat2cell(sample_table_Delta.d_gen,step_freqs_Delta); %cell array of discrete-time samples from each MCMC step

% Convert discrete-time samples into matrices of discretised probability
% distributions from each MCMC step for both variants

d_g = (0:85)';
edges = [d_g-0.5;d_g(end)+0.5];
binwidths = diff(edges);

if (max(cell2mat(days_cell_Alpha))>max(edges))||(max(cell2mat(days_cell_Delta))>max(edges))
    warning('Some data out of range')
end

counts_cell_Alpha = cellfun(@(x)histcounts(x,edges),days_cell_Alpha,'UniformOutput',false);
counts_cell_Delta = cellfun(@(x)histcounts(x,edges),days_cell_Delta,'UniformOutput',false);

counts_mat_Alpha = cell2mat(counts_cell_Alpha)';
counts_mat_Delta = cell2mat(counts_cell_Delta)';

prob_mat_Alpha = counts_mat_Alpha./(binwidths.*step_freqs_Alpha');
prob_mat_Delta = counts_mat_Delta./(binwidths.*step_freqs_Delta');

% Vectors of overall discretised probabilities for each variant

overall_probs_Alpha = mean(prob_mat_Alpha,2);
overall_probs_Delta = mean(prob_mat_Delta,2);

% Matrix of probabilities from 100 example steps of MCMC chain

steps_plot = (1:(no_steps/100):no_steps)';
prob_mat_Alpha_plot = prob_mat_Alpha(:,steps_plot);
prob_mat_Delta_plot = prob_mat_Delta(:,steps_plot);

% Plot a subset

figure(); hold on;
plot(d_g,prob_mat_Alpha_plot,'r','linewidth',0.1)
plot(d_g,overall_probs_Alpha,'k','linewidth',2)

figure(); hold on;
plot(d_g,prob_mat_Delta_plot,'r','linewidth',0.1)
plot(d_g,overall_probs_Delta,'k','linewidth',2)

% Use ksdensity to obtain line plots of continuous-time distributions

times_all_Alpha_trunc = times_all_Alpha(times_all_Alpha<=20);
times_all_Delta_trunc = times_all_Delta(times_all_Delta<=20);

t_g = 0:0.01:30;
[f_gA,~] = ksdensity(times_all_Alpha_trunc,t_g,'Support','positive','Bandwidth',0.1);
[f_gD,~] = ksdensity(times_all_Delta_trunc,t_g,'Support','positive','Bandwidth',0.1);

figure(); hold on; histogram(times_all_Alpha_trunc,'normalization','pdf')
plot(t_g,f_gA,'r','linewidth',3);
figure(); hold on; histogram(times_all_Delta_trunc,'normalization','pdf')
plot(t_g,f_gD,'r','linewidth',3);

% Save

discr_distribs.d_g = d_g;
discr_distribs.prob_mat_Alpha_plot = prob_mat_Alpha_plot;
discr_distribs.overall_probs_Alpha = overall_probs_Alpha;
discr_distribs.prob_mat_Delta_plot = prob_mat_Delta_plot;
discr_distribs.overall_probs_Delta = overall_probs_Delta;

variant_distrib_comparison.t_g = t_g;
variant_distrib_comparison.f_gA = f_gA;
variant_distrib_comparison.f_gD = f_gD;

save('../../Results/distributional_posteriors_mech.mat','discr_distribs','variant_distrib_comparison')