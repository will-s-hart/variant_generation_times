% Import posterior parameter distributions for the Alpha (denoted 1) and
% Delta (2) variants obtained by fitting the mechanistic model to household
% transmission data using data augmentation MCMC, and calculate the
% posterior distributions of the mean and standard deviation of generation
% times

% Need to run "Scripts/Parameter fitting/param_fit_mech.m" before this
% script.

clear all; close all; clc;

addpath('../../Data')
addpath('../../Results')
addpath('../../Functions/Mech')

% Load assumed parameters

load('../../Data/assumed_parameters.mat','inc_mu','inc_sigma')
load('../../Data/assumed_parameters.mat','k_inc','gamma','k_I','params_known')

% Load output of MCMC fitting procedure

load('../../Results/param_fit_mech.mat','theta_mat','ll_vec','output')

theta_mat_trace = output.theta_mat_trace; %includes burn-in, which is already removed from theta_mat

% Calculate posterior distributions of individual model parameters for each
% variant

p_E_post1 = (theta_mat(:,1));
k_E_post1 = k_inc*p_E_post1;
k_P_post1 = k_inc - k_E_post1;
mu_inv_post1 = (theta_mat(:,2));
alpha_post1 = (theta_mat(:,3));
beta0_post1 = (theta_mat(:,4));

p_E_post2 = (theta_mat(:,5));
k_E_post2 = k_inc*p_E_post2;
k_P_post2 = k_inc - k_E_post2;
mu_inv_post2 = (theta_mat(:,6));
alpha_post2 = (theta_mat(:,7));
beta0_post2 = (theta_mat(:,8));

% Posterior estimates for the mean and standard deviation of generation
% times

no_steps_kept = length(k_E_post1);
k_inc_post = repmat(params_known(1),no_steps_kept,1);
gamma_post = repmat(params_known(2),no_steps_kept,1);
k_I_post = repmat(params_known(3),no_steps_kept,1);

mu_post1 = 1./mu_inv_post1;
mu_post2 = 1./mu_inv_post2;

params_post1 = [gamma_post,mu_post1,k_inc_post,k_E_post1,k_I_post,alpha_post1,beta0_post1];
[mean_post1,sd_post1] = get_gen_mean_sd_mech(params_post1);

params_post2 = [gamma_post,mu_post2,k_inc_post,k_E_post2,k_I_post,alpha_post2,beta0_post2];
[mean_post2,sd_post2] = get_gen_mean_sd_mech(params_post2);

% Plot

figure(1); hold on;
v = violinplot([p_E_post1,p_E_post2],{'',''},'ShowData',false,'ViolinAlpha',1);
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
% ylim([0,10])
ylabel({'{\itk}_{E}/{\itk}_{inc}'})
l = legend([v1,v2],{'Alpha','Delta'},'location','n');

figure(2); hold on;
v = violinplot([mu_inv_post1,mu_inv_post2],{'',''},'ShowData',false,'ViolinAlpha',1);
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
% ylim([0,10])
ylabel({'1/\mu (days)'})
l = legend([v1,v2],{'Alpha','Delta'},'location','n');

figure(3); hold on;
v = violinplot([alpha_post1,alpha_post2],{'',''},'ShowData',false,'ViolinAlpha',1);
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
% ylim([0,10])
ylabel({'\alpha_P'})
l = legend([v1,v2],{'Alpha','Delta'},'location','n');

figure(4); hold on;
v = violinplot([beta0_post1,beta0_post2],{'',''},'ShowData',false,'ViolinAlpha',1);
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
% ylim([0,10])
ylabel({'\beta_0'})
l = legend([v1,v2],{'Alpha','Delta'},'location','n');

figure(5); hold on;
v = violinplot([mean_post1,mean_post2],{'',''},'ShowData',false,'ViolinAlpha',1);
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
% ylim([0,10])
ylabel({'Mean generation time (days)'})
l = legend([v1,v2],{'Alpha','Delta'},'location','n');

figure(6); hold on;
inds = (sd_post1<=15&sd_post2<=15); %remove outliers to ensure smooth violin
v = violinplot([sd_post1(inds),sd_post2(inds)],{'',''},'ShowData',false,'ViolinAlpha',1);
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
% ylim([0,10])
ylabel({'Standard deviation of generation';'times (days)'})
l = legend([v1,v2],{'Alpha','Delta'},'location','n');

% Save results

save('../../Results/mcmc_posterior_mech.mat','mean_post1','mean_post2','sd_post1','sd_post2','beta0_post1','beta0_post2')
save('../../Results/mcmc_posterior_mech.mat','p_E_post1','p_E_post2','mu_inv_post1','mu_inv_post2','alpha_post1','alpha_post2','-append')
save('../../Results/mcmc_posterior_mech.mat','theta_mat_trace','-append')

rmpath('../../Data')
rmpath('../../Results')
rmpath('../../Functions/Mech')