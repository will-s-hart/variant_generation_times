% Produce the panels in Figure 2 of our manuscript. This script requires
% the Violinplot-Matlab package to produce the violin plots (freely
% available at https://github.com/bastibe/Violinplot-Matlab), and the
% export-fig package to export the figure panels to PDF (freely available
% at https://github.com/altmany/export_fig)

clear all; close all; clc;

addpath('../Data')
addpath('../Results')

% Load results

load('../Results/household_gen_comp_mech2.mat','vaccine_infector_variant_comparison','vaccine_infectee_variant_comparison','vaccine_combo_variant_comparison','age_infector_variant_comparison','age_infectee_variant_comparison','date_comparison')

% Colorscheme

c1 = [0, 0.4470, 0.7410]; c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; c4 = [0.4940, 0.1840, 0.5560];
c5 = [0.4660, 0.6740, 0.1880]; c6 = [0.3010, 0.7450, 0.9330];
c7 = [0.6350, 0.0780, 0.1840]; %c8 = [0, 0.4470, 0.7410];
c8 = [1, 0.4, 0.6];

% Pre-format figures

for k = 1:6
figsetup(k)
end

% Infector vaccination status and variant

figure(1); hold on;
v = violinplot(vaccine_infector_variant_comparison,{},'ShowData',false,'ViolinAlpha',1,'BoxColor',[0,0,0]);
v(1).ViolinColor = c1; v(2).ViolinColor = c2; v(3).ViolinColor = c3; v(4).ViolinColor = c4; v(5).ViolinColor = c5; v(6).ViolinColor = c6;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot; v3 = v(3).ViolinPlot; v4 = v(4).ViolinPlot; v5 = v(5).ViolinPlot; v6 = v(6).ViolinPlot;
xticklabels([])
xlim([0.5,6.5])
ylim([0,10])
ylabel({'Mean household generation time (days)'})
l = legend([v1,v2,v3,v4,v5,v6],{'Unvaccinated, Alpha','1 dose, Alpha','2 doses, Alpha','Unvaccinated, Delta','1 dose, Delta','2 doses, Delta'});
l.FontSize = 15;
l.Position = [0.4820    0.7020    0.3600    0.2290];

% Infector vaccination status and variant

figure(2); hold on;
v = violinplot(vaccine_infectee_variant_comparison,{},'ShowData',false,'ViolinAlpha',1,'BoxColor',[0,0,0]);
v(1).ViolinColor = c1; v(2).ViolinColor = c2; v(3).ViolinColor = c3; v(4).ViolinColor = c4; v(5).ViolinColor = c5; v(6).ViolinColor = c6;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot; v3 = v(3).ViolinPlot; v4 = v(4).ViolinPlot; v5 = v(5).ViolinPlot; v6 = v(6).ViolinPlot;
xticklabels([])
xlim([0.5,6.5])
ylim([0,10])
ylabel({'Mean household generation time (days)'})
l = legend([v1,v2,v3,v4,v5,v6],{'Unvaccinated, Alpha','1 dose, Alpha','2 doses, Alpha','Unvaccinated, Delta','1 dose, Delta','2 doses, Delta'});
l.FontSize = 15;
l.Position = [0.4820    0.7020    0.3600    0.2290];

% Infector/infectee vaccination status combination and variant

figure(3); hold on;
v = violinplot(vaccine_combo_variant_comparison,{},'ShowData',false,'ViolinAlpha',1,'BoxColor',[0,0,0]);
v(1).ViolinColor = c1; v(2).ViolinColor = c2; v(3).ViolinColor = c3; v(4).ViolinColor = c4; v(5).ViolinColor = c5; v(6).ViolinColor = c6; v(7).ViolinColor = c7; v(8).ViolinColor = c8;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot; v3 = v(3).ViolinPlot; v4 = v(4).ViolinPlot; v5 = v(5).ViolinPlot; v6 = v(6).ViolinPlot; v7 = v(7).ViolinPlot; v8 = v(8).ViolinPlot;
xticklabels([])
xlim([0.5,8.5])
ylim([0,7])
ylabel({'Mean household generation time (days)'})
l = legend([v1,v2,v3,v4,v5,v6,v7,v8],{'U->U, Alpha','U->V, Alpha','V->U, Alpha','V->V, Alpha','U->U, Delta','U->V, Delta','V->U, Delta','V->V, Delta'},'NumColumns',2);
l.FontSize = 15;
l.Position = [0.3535    0.7770    0.4910    0.1550];

% Infector age and variant

figure(4); hold on;
v = violinplot(age_infector_variant_comparison,{},'ShowData',false,'ViolinAlpha',1,'BoxColor',[0,0,0]);
v(1).ViolinColor = c1; v(2).ViolinColor = c2; v(3).ViolinColor = c3; v(4).ViolinColor = c4; v(5).ViolinColor = c5; v(6).ViolinColor = c6; v(7).ViolinColor = c7; v(8).ViolinColor = c8;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot; v3 = v(3).ViolinPlot; v4 = v(4).ViolinPlot; v5 = v(5).ViolinPlot; v6 = v(6).ViolinPlot; v7 = v(7).ViolinPlot; v8 = v(8).ViolinPlot;
xticklabels([])
xlim([0.5,8.5])
ylim([0,10])
ylabel({'Mean household generation time (days)'})
l = legend([v1,v2,v3,v4,v5,v6,v7,v8],{'0-10, Alpha','11-18, Alpha','19-54, Alpha','55+, Alpha','0-10, Delta','11-18, Delta','19-54, Delta','55+, Delta'});
l.FontSize = 15;
l.Position = [0.5880    0.6110    0.2540    0.3030];

% Infectee age and variant

figure(5); hold on;
v = violinplot(age_infectee_variant_comparison,{},'ShowData',false,'ViolinAlpha',1,'BoxColor',[0,0,0]);
v(1).ViolinColor = c1; v(2).ViolinColor = c2; v(3).ViolinColor = c3; v(4).ViolinColor = c4; v(5).ViolinColor = c5; v(6).ViolinColor = c6; v(7).ViolinColor = c7; v(8).ViolinColor = c8;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot; v3 = v(3).ViolinPlot; v4 = v(4).ViolinPlot; v5 = v(5).ViolinPlot; v6 = v(6).ViolinPlot; v7 = v(7).ViolinPlot; v8 = v(8).ViolinPlot;
xticklabels([])
xlim([0.5,8.5])
ylim([0,10])
ylabel({'Mean household generation time (days)'})
l = legend([v1,v2,v3,v4,v5,v6,v7,v8],{'0-10, Alpha','11-18, Alpha','19-54, Alpha','55+, Alpha','0-10, Delta','11-18, Delta','19-54, Delta','55+, Delta'});
l.FontSize = 15;
l.Position = [0.5880    0.6110    0.2540    0.3030];

% Date

figure(6); hold on;
v = violinplot(date_comparison,{},'ShowData',false,'ViolinAlpha',1,'BoxColor',[0,0,0]);
v(1).ViolinColor = c1; v(2).ViolinColor = c2; v(3).ViolinColor = c3; v(4).ViolinColor = c4; v(5).ViolinColor = c5;  v(6).ViolinColor = c6; v(7).ViolinColor = c7;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot; v3 = v(3).ViolinPlot; v4 = v(4).ViolinPlot; v5 = v(5).ViolinPlot; v6 = v(6).ViolinPlot; v7 = v(7).ViolinPlot;
xticklabels({'Feb','Mar','Apr','May','Jun','Jul','Aug'})
xlim([0.5,7.5])
ylim([0,7])
ylabel({'Mean household generation time (days)'})

% Post-format figures

for k = 1:6
figsetup(k)
end

% Export figure panels to pdf

figure(1); export_fig Figures/Figure_2/A.pdf -nocrop -painters -transparent
figure(2); export_fig Figures/Figure_2/B.pdf -nocrop -painters -transparent
figure(3); export_fig Figures/Figure_2/C.pdf -nocrop -painters -transparent
figure(4); export_fig Figures/Figure_2/D.pdf -nocrop -painters -transparent
figure(5); export_fig Figures/Figure_2/E.pdf -nocrop -painters -transparent
figure(6); export_fig Figures/Figure_2/F.pdf -nocrop -painters -transparent

rmpath('../Data')
rmpath('../Results')