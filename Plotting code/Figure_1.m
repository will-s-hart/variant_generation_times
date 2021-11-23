% Produce the panels in Figure 1 of our manuscript. This script requires
% the Violinplot-Matlab package to produce the violin plots (freely
% available at https://github.com/bastibe/Violinplot-Matlab), and the
% export-fig package to export the figure panels to PDF (freely available
% at https://github.com/altmany/export_fig)

clear all; close all; clc;

addpath('../Data')
addpath('../Results')

% Load results

load('../Results/mcmc_posterior_mech.mat','mean_post1','mean_post2','beta0_post1','beta0_post2')
load('../Results/household_gen_comp_mech1.mat','variant_comparison','variant_sd_comparison')

% Colorscheme

c1 = [0, 0.4470, 0.7410]; c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; c4 = [0.4940, 0.1840, 0.5560];

% Pre-format figures

for k = 1:4
figsetup(k)
end

% Mean intrinsic generation time

figure(1); hold on;
v = violinplot([mean_post1,mean_post2],{'',''},'ShowData',false,'ViolinAlpha',1,'BoxColor',[0,0,0]);
v(1).ViolinColor = c1; v(2).ViolinColor = c2;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
xlim([0.5,2.5])
ylim([0,7])
ylabel({'Mean intrinsic generation time (days)'})
l = legend([v1,v2],{'Alpha','Delta'});
l.FontSize = 15;
l.Position = [0.6390    0.8090    0.1590    0.0810];

% Overall infectiousness parameter, beta_0

figure(2); hold on;
v = violinplot([beta0_post1,beta0_post2],{'',''},'ShowData',false,'ViolinAlpha',1,'BoxColor',[0,0,0]);
v(1).ViolinColor = c1; v(2).ViolinColor = c2;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
xlim([0.5,2.5])
ylim([0,6])
ylabel({'Overall transmissibility, \beta_0'})

% Mean household generation time

figure(3); hold on;
v = violinplot(variant_comparison,{},'ShowData',false,'ViolinAlpha',1,'BoxColor',[0,0,0]);
v(1).ViolinColor = c1; v(2).ViolinColor = c2;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
xlim([0.5,2.5])
xticklabels([])
ylim([0,6])
ylabel({'Mean household generation time (days)'})

% Standard deviation of household generation times

figure(4); hold on;
v = violinplot(variant_sd_comparison,{},'ShowData',false,'ViolinAlpha',1,'BoxColor',[0,0,0]);
v(1).ViolinColor = c1; v(2).ViolinColor = c2;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
xlim([0.5,2.5])
xticklabels([])
ylim([0,6])
ylabel({'Standard deviation of household';'generation times (days)'})

% Post-format figures

for k = 1:4
figsetup(k)
end

% Export figure panels to pdf

figure(1); export_fig Figures/Figure_1/A.pdf -nocrop -painters -transparent
figure(2); export_fig Figures/Figure_1/B.pdf -nocrop -painters -transparent
figure(3); export_fig Figures/Figure_1/C.pdf -nocrop -painters -transparent
figure(4); export_fig Figures/Figure_1/D.pdf -nocrop -painters -transparent

rmpath('../Data')
rmpath('../Results')