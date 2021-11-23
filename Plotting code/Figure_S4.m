% Produce the panels in Figure S4 of our manuscript. This script requires
% the export-fig package to export the figure panels to PDF (freely
% available at https://github.com/altmany/export_fig)

clear all; close all; clc;

addpath('../Data')
addpath('../Results')

% Load results

load('../Results/distributional_posteriors_mech.mat','discr_distribs','variant_distrib_comparison')

% Colorscheme

c1 = [0, 0.4470, 0.7410]; c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; c4 = [0.4940, 0.1840, 0.5560];
c5 = [0.4660, 0.6740, 0.1880]; c6 = [0.3010, 0.7450, 0.9330];
c7 = [0.6350, 0.0780, 0.1840]; %c8 = [0, 0.4470, 0.7410];
c8 = [1, 0.4, 0.6];

% Pre-format figures

for k = 1:3
figsetup(k)
end

% Alpha discretised

d_g = discr_distribs.d_g;
prob_mat_Alpha_plot = discr_distribs.prob_mat_Alpha_plot;
overall_probs_Alpha = discr_distribs.overall_probs_Alpha;

figure(1); hold on;
p1 = plot(d_g,prob_mat_Alpha_plot,'color',c6,'linewidth',0.1);
p1l = plot([-2,-1],[-2,-1],'color',c6,'linewidth',1.5); %out of sight, for legend only, to ensure visibility
p2 = plot(d_g,overall_probs_Alpha,'color',c1,'linewidth',3);
xlim([0,15])
ylim([0,0.35])
yticks(0:0.05:0.35)
xlabel({'Household generation time (days)'})
ylabel('Probability')
l = legend([p1l,p2],{'Posterior samples','Overall distribution'});
l.FontSize = 15;

% Delta discretised

prob_mat_Delta_plot = discr_distribs.prob_mat_Delta_plot;
overall_probs_Delta = discr_distribs.overall_probs_Delta;

figure(2); hold on;
p1 = plot(d_g,prob_mat_Delta_plot,'color',c8,'linewidth',0.1);
p1l = plot([-2,-1],[-2,-1],'color',c8,'linewidth',1.5); %out of sight, for legend only, to ensure visibility
p2 = plot(d_g,overall_probs_Delta,'color',c7,'linewidth',3);
xlim([0,15])
ylim([0,0.35])
yticks(0:0.05:0.35)
xlabel({'Household generation time (days)'})
ylabel('Probability')
l = legend([p1l,p2],{'Posterior samples','Overall distribution'});
l.FontSize = 15;

% Continuous time comparison

t_g = variant_distrib_comparison.t_g;
f_gA = variant_distrib_comparison.f_gA;
f_gD = variant_distrib_comparison.f_gD;

figure(3); hold on;
p1 = plot(t_g,f_gA,'color',c1,'linewidth',3);
p2 = plot(t_g,f_gD,'color',c7,'linewidth',3);
xlim([0,15])
xlabel({'Household generation time (days)'})
ylabel('Probability density')
l = legend([p1,p2],{'Alpha','Delta'});
l.FontSize = 15;

% Post-format figures

for k = 1:3
figsetup(k)
end

% Export figure panels to pdf

figure(1); export_fig Figures/Figure_S4/A.pdf -nocrop -painters -transparent
figure(2); export_fig Figures/Figure_S4/B.pdf -nocrop -painters -transparent
figure(3); export_fig Figures/Figure_S4/C.pdf -nocrop -painters -transparent

rmpath('../Data')
rmpath('../Results')