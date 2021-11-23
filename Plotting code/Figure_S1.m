% Produce the panels in Figure S4 of our manuscript. This script requires
% the export-fig package to export the figure panels to PDF (freely
% available at https://github.com/altmany/export_fig)

clear all; close all; clc;

addpath('../Data')
addpath('../Results')

% Load results

load('../Results/mcmc_posterior_mech.mat','p_E_post1','mu_inv_post1','alpha_post1','beta0_post1','p_E_post2','mu_inv_post2','alpha_post2','beta0_post2','theta_mat_trace')

% Colorscheme

c1 = [0, 0.4470, 0.7410]; c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; c4 = [0.4940, 0.1840, 0.5560];

% Pre-format figures

for k = 1:16
figsetup(k)
end

% Trace plots

theta_mat_trace = theta_mat_trace(1:10:end,:); %only keep every 1000 steps for plots (note theta_mat_trace already thinned to include only every 100 steps before this line)
no_steps_vec = 1 + 1000*((1:size(theta_mat_trace,1))'-1);

figure(1); hold on;
plot(no_steps_vec/1000000,theta_mat_trace(:,1),'color',c1,'linewidth',0.1)
ylim([0,1])
xticks(0:2.5:10)
xlabel('Step number (million)')
ylabel({'{\itk}_{E}/{\itk}_{inc}, Alpha variant'})

figure(2); hold on;
plot(no_steps_vec/1000000,theta_mat_trace(:,2),'color',c1,'linewidth',0.1)
ylim([0,10])
xticks(0:2.5:10)
xlabel('Step number (million)')
ylabel({'1/\mu (days), Alpha variant'})

figure(3); hold on;
plot(no_steps_vec/1000000,theta_mat_trace(:,3),'color',c1,'linewidth',0.1)
ylim([0,12])
xticks(0:2.5:10)
xlabel('Step number (million)')
ylabel({'\alpha_P, Alpha variant'})

figure(4); hold on;
plot(no_steps_vec/1000000,theta_mat_trace(:,4),'color',c1,'linewidth',0.1)
ylim([0,7])
xticks(0:2.5:10)
xlabel('Step number (million)')
ylabel({'\beta_0, Alpha variant'})

figure(5); hold on;
plot(no_steps_vec/1000000,theta_mat_trace(:,5),'color',c2,'linewidth',0.1)
ylim([0,1])
xticks(0:2.5:10)
xlabel('Step number (million)')
ylabel({'{\itk}_{E}/{\itk}_{inc}, Delta variant'})

figure(6); hold on;
plot(no_steps_vec/1000000,theta_mat_trace(:,6),'color',c2,'linewidth',0.1)
ylim([0,10])
xticks(0:2.5:10)
xlabel('Step number (million)')
ylabel({'1/\mu (days), Delta variant'})

figure(7); hold on;
plot(no_steps_vec/1000000,theta_mat_trace(:,7),'color',c2,'linewidth',0.1)
ylim([0,12])
xticks(0:2.5:10)
xlabel('Step number (million)')
ylabel({'\alpha_P, Delta variant'})

figure(8); hold on;
plot(no_steps_vec/1000000,theta_mat_trace(:,8),'color',c2,'linewidth',0.1)
ylim([0,7])
xticks(0:2.5:10)
xlabel('Step number (million)')
ylabel({'\beta_0, Delta variant'})

% Comparison between prior and posterior distributions

figure(9); hold on;
p2 = histogram(p_E_post1,'facecolor',c1,'normalization','pdf','edgecolor',c1);
xlim([0,1])
xticks(0:0.25:1)
p1 = plot(0:0.001:1,betapdf(0:0.001:1,2.1,2.1),'k--','linewidth',3);
ylim([0,4.5])
yticks(0:0.5:4.5)
xlabel({'{\itk}_{E}/{\itk}_{inc}, Alpha variant'})
ylabel('Density')
l = legend([p1,p2],{'Prior','Posterior'});
l.FontSize = 15;

figure(10); hold on;
p2 = histogram(mu_inv_post1,'facecolor',c1,'normalization','pdf','edgecolor',c1);
xlim([0,7])
p1 = plot(0:0.01:7,gampdf(0:0.01:7,7,0.7),'k--','linewidth',3);
xticks(1:7)
xlabel({'1/\mu (days), Alpha variant'})
ylabel('Density')

figure(11); hold on;
p2 = histogram(alpha_post1,'facecolor',c1,'normalization','pdf','edgecolor',c1);
xlim([0,10])
p1 = plot(0:0.01:10,gampdf(0:0.01:10,2.65,0.75),'k--','linewidth',3);
xticks(0:2.5:10)
xlabel({'\alpha_P, Alpha variant'})
ylabel('Density')

figure(12); hold on;
p2 = histogram(beta0_post1,'facecolor',c1,'normalization','pdf','edgecolor',c1);
xlim([0,7])
p1 = plot(0:0.01:7,gampdf(0:0.01:7,2.65,0.75),'k--','linewidth',3);
xticks(1:7)
xlabel({'\beta_0, Alpha variant'})
ylabel('Density')

figure(13); hold on;
p2 = histogram(p_E_post2,'facecolor',c2,'normalization','pdf','edgecolor',c2);
xlim([0,1])
xticks(0:0.25:1)
p1 = plot(0:0.001:1,betapdf(0:0.001:1,2.1,2.1),'k--','linewidth',3);
ylim([0,4.5])
yticks(0:0.5:4.5)
xlabel({'{\itk}_{E}/{\itk}_{inc}, Delta variant'})
ylabel('Density')
l = legend([p1,p2],{'Prior','Posterior'});
l.FontSize = 15;

figure(14); hold on;
p2 = histogram(mu_inv_post2,'facecolor',c2,'normalization','pdf','edgecolor',c2);
xlim([0,7])
p1 = plot(0:0.01:7,gampdf(0:0.01:7,7,0.7),'k--','linewidth',3);
xticks(1:7)
xlabel({'1/\mu (days), Delta variant'})
ylabel('Density')

figure(15); hold on;
p2 = histogram(alpha_post2,'facecolor',c2,'normalization','pdf','edgecolor',c2);
xlim([0,10])
p1 = plot(0:0.01:10,gampdf(0:0.01:10,2.65,0.75),'k--','linewidth',3);
xticks(0:2.5:10)
xlabel({'\alpha_P, Delta variant'})
ylabel('Density')

figure(16); hold on;
p2 = histogram(beta0_post2,'facecolor',c2,'normalization','pdf','edgecolor',c2);
xlim([0,7])
p1 = plot(0:0.01:7,gampdf(0:0.01:7,2.65,0.75),'k--','linewidth',3);
xticks(1:7)
ylim([0,2])
xlabel({'\beta_0, Delta variant'})
ylabel('Density')

% Post-format figures

for k = 1:16
figsetup(k)
end

% Export figure panels to pdf

figure(1); export_fig Figures/Figure_S1/A.pdf -nocrop -painters -transparent
figure(2); export_fig Figures/Figure_S1/B.pdf -nocrop -painters -transparent
figure(3); export_fig Figures/Figure_S1/C.pdf -nocrop -painters -transparent
figure(4); export_fig Figures/Figure_S1/D.pdf -nocrop -painters -transparent
figure(5); export_fig Figures/Figure_S1/E.pdf -nocrop -painters -transparent
figure(6); export_fig Figures/Figure_S1/F.pdf -nocrop -painters -transparent
figure(7); export_fig Figures/Figure_S1/G.pdf -nocrop -painters -transparent
figure(8); export_fig Figures/Figure_S1/H.pdf -nocrop -painters -transparent

figure(9); export_fig Figures/Figure_S1/I.pdf -nocrop -painters -transparent
figure(10); export_fig Figures/Figure_S1/J.pdf -nocrop -painters -transparent
figure(11); export_fig Figures/Figure_S1/K.pdf -nocrop -painters -transparent
figure(12); export_fig Figures/Figure_S1/L.pdf -nocrop -painters -transparent
figure(13); export_fig Figures/Figure_S1/M.pdf -nocrop -painters -transparent
figure(14); export_fig Figures/Figure_S1/N.pdf -nocrop -painters -transparent
figure(15); export_fig Figures/Figure_S1/O.pdf -nocrop -painters -transparent
figure(16); export_fig Figures/Figure_S1/P.pdf -nocrop -painters -transparent

rmpath('../Data')
rmpath('../Results')