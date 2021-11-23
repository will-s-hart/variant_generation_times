% Produce the panels in Figure S3 of our manuscript. This script requires
% the Violinplot-Matlab package to produce the violin plots (freely
% available at https://github.com/bastibe/Violinplot-Matlab), and the
% export-fig package to export the figure panels to PDF (freely available
% at https://github.com/altmany/export_fig)

clear all; close all; clc;

addpath('../Data')
addpath('../Results')

% Load results

load('../Results/mcmc_posterior_mech.mat','p_E_post1','p_E_post2','mu_inv_post1','mu_inv_post2','alpha_post1','alpha_post2','sd_post1','sd_post2')

% Colorscheme

c1 = [0, 0.4470, 0.7410]; c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; c4 = [0.4940, 0.1840, 0.5560];

% Pre-format figures

for k = 1:4
figsetup(k)
end

% Ratio between the shape parameters of the Gamma distributed latent (E)
% and incubation (combined E and P) periods, p_E=k_E/k_inc

figure(1); hold on;
v = violinplot([p_E_post1,p_E_post2],{'',''},'ShowData',false,'ViolinAlpha',1,'BoxColor',[0,0,0]);
v(1).ViolinColor = c1; v(2).ViolinColor = c2;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
xlim([0.5,2.5])
ylim([0,1])
ylabel({'{\itk}_{E}/{\itk}_{inc}'})
l = legend([v1,v2],{'Alpha','Delta'});
l.FontSize = 15;

% Mean symptomatic infectious (I) period, mu_inv=1/mu

figure(2); hold on;
v = violinplot([mu_inv_post1,mu_inv_post2],{'',''},'ShowData',false,'ViolinAlpha',1,'BoxColor',[0,0,0]);
v(1).ViolinColor = c1; v(2).ViolinColor = c2;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
xlim([0.5,2.5])
ylim([0,7])
ylabel({'1/\mu (days)'})

% Ratio of the transmission rates during the presymptomatic infectious (P)
% and symptomatic infectious (I) periods, alpha=alpha_P

figure(3); hold on;
inds = (alpha_post1<=15&alpha_post2<=15);
v = violinplot([alpha_post1(inds),alpha_post2(inds)],{'',''},'ShowData',false,'ViolinAlpha',1,'BoxColor',[0,0,0]);
v(1).ViolinColor = c1; v(2).ViolinColor = c2;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
xlim([0.5,2.5])
ylim([0,15])
ylabel({'\alpha_P'})

% Standard deviation of intrinsic generation time distribution

figure(4); hold on;
v = violinplot([sd_post1,sd_post2],{'',''},'ShowData',false,'ViolinAlpha',1,'BoxColor',[0,0,0]);
v(1).ViolinColor = c1; v(2).ViolinColor = c2;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
xlim([0.5,2.5])
ylim([0,6])
ylabel({'Standard deviation of intrinsic';'generation time distribution (days)'})

% Post-format figures

for k = 1:4
figsetup(k)
end

% Export figure panels to pdf

figure(1); export_fig Figures/Figure_S3/A.pdf -nocrop -painters -transparent
figure(2); export_fig Figures/Figure_S3/B.pdf -nocrop -painters -transparent
figure(3); export_fig Figures/Figure_S3/C.pdf -nocrop -painters -transparent
figure(4); export_fig Figures/Figure_S3/D.pdf -nocrop -painters -transparent

rmpath('../Data')
rmpath('../Results')