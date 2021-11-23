% Compare posterior estimates of the mean household generation time by
% infector/infectee/household variant, vaccination status, age and date,
% and compare the standard deviation by variant

% Need to run "Scripts/Parameter fitting/hh_gen_samples_mech.m" before this
% script.

clear all; close all; clc;

addpath('../../Data')

% Load data

load('../../Data/data.mat','data_struct_observed')

% Load output of MCMC fitting procedure

load('../../Results/hh_gen_samples_mech.mat','sample_table','sample_table_all')

% Vector of edges used to determine number of GT samples in different
% categories per step of the MCMC chain

no_steps = size(sample_table,1);
edges = (0.5:(no_steps+1))';

% Compare based on variant

variant = data_struct_observed.variant;

Alpha = (variant==1);
Delta = (variant==2);

Alpha_inds = Alpha(sample_table_all.infector);
Delta_inds = Delta(sample_table_all.infector);

sample_table_Alpha = sample_table_all(Alpha_inds,:);
sample_table_Delta = sample_table_all(Delta_inds,:);

% Mean when compared by variant

[step_freqs_Alpha,~] = histcounts(sample_table_Alpha.step_no,edges);
step_freqs_Alpha = step_freqs_Alpha';

times_cell_Alpha = mat2cell(sample_table_Alpha.t_gen,step_freqs_Alpha);
mean_post_Alpha = cellfun(@mean,times_cell_Alpha);

[step_freqs_Delta,~] = histcounts(sample_table_Delta.step_no,edges);
step_freqs_Delta = step_freqs_Delta';

times_cell_Delta = mat2cell(sample_table_Delta.t_gen,step_freqs_Delta);
mean_post_Delta = cellfun(@mean,times_cell_Delta);

figure(); hold on;
v = violinplot([mean_post_Alpha,mean_post_Delta],{'',''},'ShowData',false,'ViolinAlpha',1);
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
ylabel({'Mean generation time (days)'})
l = legend([v1,v2],{'Alpha','Delta'});
ylim([0,6])

% SD when compared by variant

sd_post_Alpha = cellfun(@std,times_cell_Alpha);
sd_post_Delta = cellfun(@std,times_cell_Delta);

figure(); hold on;
v = violinplot([sd_post_Alpha,sd_post_Delta],{'',''},'ShowData',false,'ViolinAlpha',1);
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
ylabel({'Standard deviation of generation times (days)'})
l = legend([v1,v2],{'Alpha','Delta'});
ylim([0,5])

% Compare based on both variant and vaccine status of infector

vaccine = data_struct_observed.vaccine;
unvacc = (vaccine==0);
onedose = (vaccine==1);
twodose = (vaccine==2);

unvacc_Alpha = unvacc&Alpha;
onedose_Alpha = onedose&Alpha;
twodose_Alpha = twodose&Alpha;
unvacc_Delta = unvacc&Delta;
onedose_Delta = onedose&Delta;
twodose_Delta = twodose&Delta;

unvacc_Alpha_inds = unvacc_Alpha(sample_table_all.infector);
onedose_Alpha_inds = onedose_Alpha(sample_table_all.infector);
twodose_Alpha_inds = twodose_Alpha(sample_table_all.infector);
unvacc_Delta_inds = unvacc_Delta(sample_table_all.infector);
onedose_Delta_inds = onedose_Delta(sample_table_all.infector);
twodose_Delta_inds = twodose_Delta(sample_table_all.infector);

sample_table_unvacc_Alpha = sample_table_all(unvacc_Alpha_inds,:);
sample_table_onedose_Alpha = sample_table_all(onedose_Alpha_inds,:);
sample_table_twodose_Alpha = sample_table_all(twodose_Alpha_inds,:);
sample_table_unvacc_Delta = sample_table_all(unvacc_Delta_inds,:);
sample_table_onedose_Delta = sample_table_all(onedose_Delta_inds,:);
sample_table_twodose_Delta = sample_table_all(twodose_Delta_inds,:);

% Mean when compared by vaccine/variant status of infector

[step_freqs_unvacc_Alpha,~] = histcounts(sample_table_unvacc_Alpha.step_no,edges);
step_freqs_unvacc_Alpha = step_freqs_unvacc_Alpha';

times_cell_unvacc_Alpha = mat2cell(sample_table_unvacc_Alpha.t_gen,step_freqs_unvacc_Alpha);
mean_post_unvacc_Alpha = cellfun(@mean,times_cell_unvacc_Alpha);

[step_freqs_onedose_Alpha,~] = histcounts(sample_table_onedose_Alpha.step_no,edges);
step_freqs_onedose_Alpha = step_freqs_onedose_Alpha';

times_cell_onedose_Alpha = mat2cell(sample_table_onedose_Alpha.t_gen,step_freqs_onedose_Alpha);
mean_post_onedose_Alpha = cellfun(@mean,times_cell_onedose_Alpha);

[step_freqs_twodose_Alpha,~] = histcounts(sample_table_twodose_Alpha.step_no,edges);
step_freqs_twodose_Alpha = step_freqs_twodose_Alpha';

times_cell_twodose_Alpha = mat2cell(sample_table_twodose_Alpha.t_gen,step_freqs_twodose_Alpha);
mean_post_twodose_Alpha = cellfun(@mean,times_cell_twodose_Alpha);

[step_freqs_unvacc_Delta,~] = histcounts(sample_table_unvacc_Delta.step_no,edges);
step_freqs_unvacc_Delta = step_freqs_unvacc_Delta';

times_cell_unvacc_Delta = mat2cell(sample_table_unvacc_Delta.t_gen,step_freqs_unvacc_Delta);
mean_post_unvacc_Delta = cellfun(@mean,times_cell_unvacc_Delta);

[step_freqs_onedose_Delta,~] = histcounts(sample_table_onedose_Delta.step_no,edges);
step_freqs_onedose_Delta = step_freqs_onedose_Delta';

times_cell_onedose_Delta = mat2cell(sample_table_onedose_Delta.t_gen,step_freqs_onedose_Delta);
mean_post_onedose_Delta = cellfun(@mean,times_cell_onedose_Delta);

[step_freqs_twodose_Delta,~] = histcounts(sample_table_twodose_Delta.step_no,edges);
step_freqs_twodose_Delta = step_freqs_twodose_Delta';

times_cell_twodose_Delta = mat2cell(sample_table_twodose_Delta.t_gen,step_freqs_twodose_Delta);
mean_post_twodose_Delta = cellfun(@mean,times_cell_twodose_Delta);

figure(); hold on;
v = violinplot([mean_post_unvacc_Alpha,mean_post_onedose_Alpha,mean_post_twodose_Alpha,mean_post_unvacc_Delta,mean_post_onedose_Delta,mean_post_twodose_Delta],{'','','','','',''},'ShowData',false,'ViolinAlpha',1);
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot; v3 = v(3).ViolinPlot; v4 = v(4).ViolinPlot; v5 = v(5).ViolinPlot; v6 = v(6).ViolinPlot;
ylabel({'Mean generation time (days)'})
l = legend([v1,v2,v3,v4,v5,v6],{'0 doses, Alpha','1 dose, Alpha','2 doses, Alpha','0 doses, Delta','1 dose, Delta','2 doses, Delta'});
ylim([0,10])

% Compare based on both variant and vaccine status of infectee

unvacc_infectee_Alpha_inds = unvacc_Alpha(sample_table_all.infectee);
onedose_infectee_Alpha_inds = onedose_Alpha(sample_table_all.infectee);
twodose_infectee_Alpha_inds = twodose_Alpha(sample_table_all.infectee);
unvacc_infectee_Delta_inds = unvacc_Delta(sample_table_all.infectee);
onedose_infectee_Delta_inds = onedose_Delta(sample_table_all.infectee);
twodose_infectee_Delta_inds = twodose_Delta(sample_table_all.infectee);

sample_table_unvacc_infectee_Alpha = sample_table_all(unvacc_infectee_Alpha_inds,:);
sample_table_onedose_infectee_Alpha = sample_table_all(onedose_infectee_Alpha_inds,:);
sample_table_twodose_infectee_Alpha = sample_table_all(twodose_infectee_Alpha_inds,:);
sample_table_unvacc_infectee_Delta = sample_table_all(unvacc_infectee_Delta_inds,:);
sample_table_onedose_infectee_Delta = sample_table_all(onedose_infectee_Delta_inds,:);
sample_table_twodose_infectee_Delta = sample_table_all(twodose_infectee_Delta_inds,:);

% Mean when compared by vaccine/variant status of infectee

[step_freqs_unvacc_infectee_Alpha,~] = histcounts(sample_table_unvacc_infectee_Alpha.step_no,edges);
step_freqs_unvacc_infectee_Alpha = step_freqs_unvacc_infectee_Alpha';

times_cell_unvacc_infectee_Alpha = mat2cell(sample_table_unvacc_infectee_Alpha.t_gen,step_freqs_unvacc_infectee_Alpha);
mean_post_unvacc_infectee_Alpha = cellfun(@mean,times_cell_unvacc_infectee_Alpha);

[step_freqs_onedose_infectee_Alpha,~] = histcounts(sample_table_onedose_infectee_Alpha.step_no,edges);
step_freqs_onedose_infectee_Alpha = step_freqs_onedose_infectee_Alpha';

times_cell_onedose_infectee_Alpha = mat2cell(sample_table_onedose_infectee_Alpha.t_gen,step_freqs_onedose_infectee_Alpha);
mean_post_onedose_infectee_Alpha = cellfun(@mean,times_cell_onedose_infectee_Alpha);

[step_freqs_twodose_infectee_Alpha,~] = histcounts(sample_table_twodose_infectee_Alpha.step_no,edges);
step_freqs_twodose_infectee_Alpha = step_freqs_twodose_infectee_Alpha';

times_cell_twodose_infectee_Alpha = mat2cell(sample_table_twodose_infectee_Alpha.t_gen,step_freqs_twodose_infectee_Alpha);
mean_post_twodose_infectee_Alpha = cellfun(@mean,times_cell_twodose_infectee_Alpha);

[step_freqs_unvacc_infectee_Delta,~] = histcounts(sample_table_unvacc_infectee_Delta.step_no,edges);
step_freqs_unvacc_infectee_Delta = step_freqs_unvacc_infectee_Delta';

times_cell_unvacc_infectee_Delta = mat2cell(sample_table_unvacc_infectee_Delta.t_gen,step_freqs_unvacc_infectee_Delta);
mean_post_unvacc_infectee_Delta = cellfun(@mean,times_cell_unvacc_infectee_Delta);

[step_freqs_onedose_infectee_Delta,~] = histcounts(sample_table_onedose_infectee_Delta.step_no,edges);
step_freqs_onedose_infectee_Delta = step_freqs_onedose_infectee_Delta';

times_cell_onedose_infectee_Delta = mat2cell(sample_table_onedose_infectee_Delta.t_gen,step_freqs_onedose_infectee_Delta);
mean_post_onedose_infectee_Delta = cellfun(@mean,times_cell_onedose_infectee_Delta);

[step_freqs_twodose_infectee_Delta,~] = histcounts(sample_table_twodose_infectee_Delta.step_no,edges);
step_freqs_twodose_infectee_Delta = step_freqs_twodose_infectee_Delta';

times_cell_twodose_infectee_Delta = mat2cell(sample_table_twodose_infectee_Delta.t_gen,step_freqs_twodose_infectee_Delta);
mean_post_twodose_infectee_Delta = cellfun(@mean,times_cell_twodose_infectee_Delta);

figure(); hold on;
v = violinplot([mean_post_unvacc_infectee_Alpha,mean_post_onedose_infectee_Alpha,mean_post_twodose_infectee_Alpha,mean_post_unvacc_infectee_Delta,mean_post_onedose_infectee_Delta,mean_post_twodose_infectee_Delta],{'','','','','',''},'ShowData',false,'ViolinAlpha',1);
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot; v3 = v(3).ViolinPlot; v4 = v(4).ViolinPlot; v5 = v(5).ViolinPlot; v6 = v(6).ViolinPlot;
ylabel({'Mean generation time (days)'})
l = legend([v1,v2,v3,v4,v5,v6],{'0 doses, Alpha','1 dose, Alpha','2 doses, Alpha','0 doses, Delta','1 dose, Delta','2 doses, Delta'});
ylim([0,10])

% Compare based on both variant and vaccine status of infector/infectee

vacc_Alpha = (onedose_Alpha|twodose_Alpha);
vacc_Delta = (onedose_Delta|twodose_Delta);

uunvacc_Alpha_inds = unvacc_Alpha(sample_table_all.infector)&unvacc_Alpha(sample_table_all.infectee);
uvacc_Alpha_inds = unvacc_Alpha(sample_table_all.infector)&vacc_Alpha(sample_table_all.infectee);
vunvacc_Alpha_inds = vacc_Alpha(sample_table_all.infector)&unvacc_Alpha(sample_table_all.infectee);
vvacc_Alpha_inds = vacc_Alpha(sample_table_all.infector)&vacc_Alpha(sample_table_all.infectee);
uu_Delta_inds = unvacc_Delta(sample_table_all.infector)&unvacc_Delta(sample_table_all.infectee);
uv_Delta_inds = unvacc_Delta(sample_table_all.infector)&vacc_Delta(sample_table_all.infectee);
vu_Delta_inds = vacc_Delta(sample_table_all.infector)&unvacc_Delta(sample_table_all.infectee);
vv_Delta_inds = vacc_Delta(sample_table_all.infector)&vacc_Delta(sample_table_all.infectee);

sample_table_uunvacc_Alpha = sample_table_all(uunvacc_Alpha_inds,:);
sample_table_uvacc_Alpha = sample_table_all(uvacc_Alpha_inds,:);
sample_table_vunvacc_Alpha = sample_table_all(vunvacc_Alpha_inds,:);
sample_table_vvacc_Alpha = sample_table_all(vvacc_Alpha_inds,:);
sample_table_uu_Delta = sample_table_all(uu_Delta_inds,:);
sample_table_uv_Delta = sample_table_all(uv_Delta_inds,:);
sample_table_vu_Delta = sample_table_all(vu_Delta_inds,:);
sample_table_vv_Delta = sample_table_all(vv_Delta_inds,:);

% Mean when compared by vaccine/variant status of infector/infectee

[step_freqs_uunvacc_Alpha,~] = histcounts(sample_table_uunvacc_Alpha.step_no,edges);
step_freqs_uunvacc_Alpha = step_freqs_uunvacc_Alpha';

times_cell_uunvacc_Alpha = mat2cell(sample_table_uunvacc_Alpha.t_gen,step_freqs_uunvacc_Alpha);
mean_post_uunvacc_Alpha = cellfun(@mean,times_cell_uunvacc_Alpha);

[step_freqs_uvacc_Alpha,~] = histcounts(sample_table_uvacc_Alpha.step_no,edges);
step_freqs_uvacc_Alpha = step_freqs_uvacc_Alpha';

times_cell_uvacc_Alpha = mat2cell(sample_table_uvacc_Alpha.t_gen,step_freqs_uvacc_Alpha);
mean_post_uvacc_Alpha = cellfun(@mean,times_cell_uvacc_Alpha);

[step_freqs_vunvacc_Alpha,~] = histcounts(sample_table_vunvacc_Alpha.step_no,edges);
step_freqs_vunvacc_Alpha = step_freqs_vunvacc_Alpha';

times_cell_vunvacc_Alpha = mat2cell(sample_table_vunvacc_Alpha.t_gen,step_freqs_vunvacc_Alpha);
mean_post_vunvacc_Alpha = cellfun(@mean,times_cell_vunvacc_Alpha);

[step_freqs_vvacc_Alpha,~] = histcounts(sample_table_vvacc_Alpha.step_no,edges);
step_freqs_vvacc_Alpha = step_freqs_vvacc_Alpha';

times_cell_vvacc_Alpha = mat2cell(sample_table_vvacc_Alpha.t_gen,step_freqs_vvacc_Alpha);
mean_post_vvacc_Alpha = cellfun(@mean,times_cell_vvacc_Alpha);

[step_freqs_uu_Delta,~] = histcounts(sample_table_uu_Delta.step_no,edges);
step_freqs_uu_Delta = step_freqs_uu_Delta';

times_cell_uu_Delta = mat2cell(sample_table_uu_Delta.t_gen,step_freqs_uu_Delta);
mean_post_uu_Delta = cellfun(@mean,times_cell_uu_Delta);

[step_freqs_uv_Delta,~] = histcounts(sample_table_uv_Delta.step_no,edges);
step_freqs_uv_Delta = step_freqs_uv_Delta';

times_cell_uv_Delta = mat2cell(sample_table_uv_Delta.t_gen,step_freqs_uv_Delta);
mean_post_uv_Delta = cellfun(@mean,times_cell_uv_Delta);

[step_freqs_vu_Delta,~] = histcounts(sample_table_vu_Delta.step_no,edges);
step_freqs_vu_Delta = step_freqs_vu_Delta';

times_cell_vu_Delta = mat2cell(sample_table_vu_Delta.t_gen,step_freqs_vu_Delta);
mean_post_vu_Delta = cellfun(@mean,times_cell_vu_Delta);

[step_freqs_vv_Delta,~] = histcounts(sample_table_vv_Delta.step_no,edges);
step_freqs_vv_Delta = step_freqs_vv_Delta';

times_cell_vv_Delta = mat2cell(sample_table_vv_Delta.t_gen,step_freqs_vv_Delta);
mean_post_vv_Delta = cellfun(@mean,times_cell_vv_Delta);

figure(); hold on;
v = violinplot([mean_post_uunvacc_Alpha,mean_post_uvacc_Alpha,mean_post_vunvacc_Alpha,mean_post_vvacc_Alpha,mean_post_uu_Delta,mean_post_uv_Delta,mean_post_vu_Delta,mean_post_vv_Delta],{'','','','','','','',''},'ShowData',false,'ViolinAlpha',1);
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot; v3 = v(3).ViolinPlot; v4 = v(4).ViolinPlot; v5 = v(5).ViolinPlot; v6 = v(6).ViolinPlot; v7 = v(7).ViolinPlot; v8 = v(8).ViolinPlot;
ylabel({'Mean generation time (days)'})
l = legend([v1,v2,v3,v4,v5,v6,v7,v8],{'Unvacc->unvacc, Alpha','Unvacc->vacc, Alpha','Vacc->unvacc, Alpha','Vacc->vacc, Alpha','Unvacc->unvacc, Delta','Unvacc->vacc, Delta','Vacc->unvacc, Delta','Vacc->vacc, Delta'});
ylim([0,10])

% Compare by age

age_group = data_struct_observed.age_group;
age1 = (age_group==1);
age2 = (age_group==2);
age3 = (age_group==3);
age4 = (age_group==4);

age1_Alpha = age1&Alpha;
age2_Alpha = age2&Alpha;
age3_Alpha = age3&Alpha;
age4_Alpha = age4&Alpha;
age1_Delta = age1&Delta;
age2_Delta = age2&Delta;
age3_Delta = age3&Delta;
age4_Delta = age4&Delta;

age1_Alpha_inds = age1_Alpha(sample_table_all.infector);
age2_Alpha_inds = age2_Alpha(sample_table_all.infector);
age3_Alpha_inds = age3_Alpha(sample_table_all.infector);
age4_Alpha_inds = age4_Alpha(sample_table_all.infector);
age1_Delta_inds = age1_Delta(sample_table_all.infector);
age2_Delta_inds = age2_Delta(sample_table_all.infector);
age3_Delta_inds = age3_Delta(sample_table_all.infector);
age4_Delta_inds = age4_Delta(sample_table_all.infector);

sample_table_age1_Alpha = sample_table_all(age1_Alpha_inds,:);
sample_table_age2_Alpha = sample_table_all(age2_Alpha_inds,:);
sample_table_age3_Alpha = sample_table_all(age3_Alpha_inds,:);
sample_table_age4_Alpha = sample_table_all(age4_Alpha_inds,:);
sample_table_age1_Delta = sample_table_all(age1_Delta_inds,:);
sample_table_age2_Delta = sample_table_all(age2_Delta_inds,:);
sample_table_age3_Delta = sample_table_all(age3_Delta_inds,:);
sample_table_age4_Delta = sample_table_all(age4_Delta_inds,:);

% Mean when compared by age and variant

[step_freqs_age1_Alpha,~] = histcounts(sample_table_age1_Alpha.step_no,edges);
step_freqs_age1_Alpha = step_freqs_age1_Alpha';

times_cell_age1_Alpha = mat2cell(sample_table_age1_Alpha.t_gen,step_freqs_age1_Alpha);
mean_post_age1_Alpha = cellfun(@mean,times_cell_age1_Alpha);

[step_freqs_age2_Alpha,~] = histcounts(sample_table_age2_Alpha.step_no,edges);
step_freqs_age2_Alpha = step_freqs_age2_Alpha';

times_cell_age2_Alpha = mat2cell(sample_table_age2_Alpha.t_gen,step_freqs_age2_Alpha);
mean_post_age2_Alpha = cellfun(@mean,times_cell_age2_Alpha);

[step_freqs_age3_Alpha,~] = histcounts(sample_table_age3_Alpha.step_no,edges);
step_freqs_age3_Alpha = step_freqs_age3_Alpha';

times_cell_age3_Alpha = mat2cell(sample_table_age3_Alpha.t_gen,step_freqs_age3_Alpha);
mean_post_age3_Alpha = cellfun(@mean,times_cell_age3_Alpha);

[step_freqs_age4_Alpha,~] = histcounts(sample_table_age4_Alpha.step_no,edges);
step_freqs_age4_Alpha = step_freqs_age4_Alpha';

times_cell_age4_Alpha = mat2cell(sample_table_age4_Alpha.t_gen,step_freqs_age4_Alpha);
mean_post_age4_Alpha = cellfun(@mean,times_cell_age4_Alpha);

[step_freqs_age1_Delta,~] = histcounts(sample_table_age1_Delta.step_no,edges);
step_freqs_age1_Delta = step_freqs_age1_Delta';

times_cell_age1_Delta = mat2cell(sample_table_age1_Delta.t_gen,step_freqs_age1_Delta);
mean_post_age1_Delta = cellfun(@mean,times_cell_age1_Delta);

[step_freqs_age2_Delta,~] = histcounts(sample_table_age2_Delta.step_no,edges);
step_freqs_age2_Delta = step_freqs_age2_Delta';

times_cell_age2_Delta = mat2cell(sample_table_age2_Delta.t_gen,step_freqs_age2_Delta);
mean_post_age2_Delta = cellfun(@mean,times_cell_age2_Delta);

[step_freqs_age3_Delta,~] = histcounts(sample_table_age3_Delta.step_no,edges);
step_freqs_age3_Delta = step_freqs_age3_Delta';

times_cell_age3_Delta = mat2cell(sample_table_age3_Delta.t_gen,step_freqs_age3_Delta);
mean_post_age3_Delta = cellfun(@mean,times_cell_age3_Delta);

[step_freqs_age4_Delta,~] = histcounts(sample_table_age4_Delta.step_no,edges);
step_freqs_age4_Delta = step_freqs_age4_Delta';

times_cell_age4_Delta = mat2cell(sample_table_age4_Delta.t_gen,step_freqs_age4_Delta);
mean_post_age4_Delta = cellfun(@mean,times_cell_age4_Delta);

figure(); hold on;
v = violinplot([mean_post_age1_Alpha,mean_post_age2_Alpha,mean_post_age3_Alpha,mean_post_age4_Alpha,mean_post_age1_Delta,mean_post_age2_Delta,mean_post_age3_Delta,mean_post_age4_Delta],{'','','','','','','',''},'ShowData',false,'ViolinAlpha',1);
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot; v3 = v(3).ViolinPlot; v4 = v(4).ViolinPlot; v5 = v(5).ViolinPlot; v6 = v(6).ViolinPlot; v7 = v(7).ViolinPlot; v8 = v(8).ViolinPlot;
ylabel({'Mean generation time (days)'})
l = legend([v1,v2,v3,v4,v5,v6,v7,v8],{'0-10, Alpha','11-18, Alpha','19-54, Alpha','55+, Alpha','0-10, Delta','11-18, Delta','19-54, Delta','55+, Delta'});
ylim([0,10])

% Compare based on age of infectee and variant

age1_infectee_Alpha_inds = age1_Alpha(sample_table_all.infectee);
age2_infectee_Alpha_inds = age2_Alpha(sample_table_all.infectee);
age3_infectee_Alpha_inds = age3_Alpha(sample_table_all.infectee);
age4_infectee_Alpha_inds = age4_Alpha(sample_table_all.infectee);
age1_infectee_Delta_inds = age1_Delta(sample_table_all.infectee);
age2_infectee_Delta_inds = age2_Delta(sample_table_all.infectee);
age3_infectee_Delta_inds = age3_Delta(sample_table_all.infectee);
age4_infectee_Delta_inds = age4_Delta(sample_table_all.infectee);

sample_table_age1_infectee_Alpha = sample_table_all(age1_infectee_Alpha_inds,:);
sample_table_age2_infectee_Alpha = sample_table_all(age2_infectee_Alpha_inds,:);
sample_table_age3_infectee_Alpha = sample_table_all(age3_infectee_Alpha_inds,:);
sample_table_age4_infectee_Alpha = sample_table_all(age4_infectee_Alpha_inds,:);
sample_table_age1_infectee_Delta = sample_table_all(age1_infectee_Delta_inds,:);
sample_table_age2_infectee_Delta = sample_table_all(age2_infectee_Delta_inds,:);
sample_table_age3_infectee_Delta = sample_table_all(age3_infectee_Delta_inds,:);
sample_table_age4_infectee_Delta = sample_table_all(age4_infectee_Delta_inds,:);

% Mean when compared by age of infectee and variant

[step_freqs_age1_infectee_Alpha,~] = histcounts(sample_table_age1_infectee_Alpha.step_no,edges);
step_freqs_age1_infectee_Alpha = step_freqs_age1_infectee_Alpha';

times_cell_age1_infectee_Alpha = mat2cell(sample_table_age1_infectee_Alpha.t_gen,step_freqs_age1_infectee_Alpha);
mean_post_age1_infectee_Alpha = cellfun(@mean,times_cell_age1_infectee_Alpha);

[step_freqs_age2_infectee_Alpha,~] = histcounts(sample_table_age2_infectee_Alpha.step_no,edges);
step_freqs_age2_infectee_Alpha = step_freqs_age2_infectee_Alpha';

times_cell_age2_infectee_Alpha = mat2cell(sample_table_age2_infectee_Alpha.t_gen,step_freqs_age2_infectee_Alpha);
mean_post_age2_infectee_Alpha = cellfun(@mean,times_cell_age2_infectee_Alpha);

[step_freqs_age3_infectee_Alpha,~] = histcounts(sample_table_age3_infectee_Alpha.step_no,edges);
step_freqs_age3_infectee_Alpha = step_freqs_age3_infectee_Alpha';

times_cell_age3_infectee_Alpha = mat2cell(sample_table_age3_infectee_Alpha.t_gen,step_freqs_age3_infectee_Alpha);
mean_post_age3_infectee_Alpha = cellfun(@mean,times_cell_age3_infectee_Alpha);

[step_freqs_age4_infectee_Alpha,~] = histcounts(sample_table_age4_infectee_Alpha.step_no,edges);
step_freqs_age4_infectee_Alpha = step_freqs_age4_infectee_Alpha';

times_cell_age4_infectee_Alpha = mat2cell(sample_table_age4_infectee_Alpha.t_gen,step_freqs_age4_infectee_Alpha);
mean_post_age4_infectee_Alpha = cellfun(@mean,times_cell_age4_infectee_Alpha);

[step_freqs_age1_infectee_Delta,~] = histcounts(sample_table_age1_infectee_Delta.step_no,edges);
step_freqs_age1_infectee_Delta = step_freqs_age1_infectee_Delta';

times_cell_age1_infectee_Delta = mat2cell(sample_table_age1_infectee_Delta.t_gen,step_freqs_age1_infectee_Delta);
mean_post_age1_infectee_Delta = cellfun(@mean,times_cell_age1_infectee_Delta);

[step_freqs_age2_infectee_Delta,~] = histcounts(sample_table_age2_infectee_Delta.step_no,edges);
step_freqs_age2_infectee_Delta = step_freqs_age2_infectee_Delta';

times_cell_age2_infectee_Delta = mat2cell(sample_table_age2_infectee_Delta.t_gen,step_freqs_age2_infectee_Delta);
mean_post_age2_infectee_Delta = cellfun(@mean,times_cell_age2_infectee_Delta);

[step_freqs_age3_infectee_Delta,~] = histcounts(sample_table_age3_infectee_Delta.step_no,edges);
step_freqs_age3_infectee_Delta = step_freqs_age3_infectee_Delta';

times_cell_age3_infectee_Delta = mat2cell(sample_table_age3_infectee_Delta.t_gen,step_freqs_age3_infectee_Delta);
mean_post_age3_infectee_Delta = cellfun(@mean,times_cell_age3_infectee_Delta);

[step_freqs_age4_infectee_Delta,~] = histcounts(sample_table_age4_infectee_Delta.step_no,edges);
step_freqs_age4_infectee_Delta = step_freqs_age4_infectee_Delta';

times_cell_age4_infectee_Delta = mat2cell(sample_table_age4_infectee_Delta.t_gen,step_freqs_age4_infectee_Delta);
mean_post_age4_infectee_Delta = cellfun(@mean,times_cell_age4_infectee_Delta);

figure(); hold on;
v = violinplot([mean_post_age1_infectee_Alpha,mean_post_age2_infectee_Alpha,mean_post_age3_infectee_Alpha,mean_post_age4_infectee_Alpha,mean_post_age1_infectee_Delta,mean_post_age2_infectee_Delta,mean_post_age3_infectee_Delta,mean_post_age4_infectee_Delta],{'','','','','','','',''},'ShowData',false,'ViolinAlpha',1);
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot; v3 = v(3).ViolinPlot; v4 = v(4).ViolinPlot; v5 = v(5).ViolinPlot; v6 = v(6).ViolinPlot; v7 = v(7).ViolinPlot; v8 = v(8).ViolinPlot;
ylabel({'Mean generation time (days)'})
l = legend([v1,v2,v3,v4,v5,v6,v7,v8],{'0-10, Alpha','11-18, Alpha','19-54, Alpha','55+, Alpha','0-10, Delta','11-18, Delta','19-54, Delta','55+, Delta'});
ylim([0,10])

% Compare by month

recruitment_month = data_struct_observed.recruitment_month;

feb = (recruitment_month==2);
mar = (recruitment_month==3);
apr = (recruitment_month==4);
may = (recruitment_month==5);
jun = (recruitment_month==6);
jul = (recruitment_month==7);
aug = (recruitment_month==8);

feb_inds = feb(sample_table_all.infector);
mar_inds = mar(sample_table_all.infector);
apr_inds = apr(sample_table_all.infector);
may_inds = may(sample_table_all.infector);
jun_inds = jun(sample_table_all.infector);
jul_inds = jul(sample_table_all.infector);
aug_inds = aug(sample_table_all.infector);

sample_table_feb = sample_table_all(feb_inds,:);
sample_table_mar = sample_table_all(mar_inds,:);
sample_table_apr = sample_table_all(apr_inds,:);
sample_table_may = sample_table_all(may_inds,:);
sample_table_jun = sample_table_all(jun_inds,:);
sample_table_jul = sample_table_all(jul_inds,:);
sample_table_aug = sample_table_all(aug_inds,:);

% Mean when compared by month

[step_freqs_feb,~] = histcounts(sample_table_feb.step_no,edges);
step_freqs_feb = step_freqs_feb';

times_cell_feb = mat2cell(sample_table_feb.t_gen,step_freqs_feb);
mean_post_feb = cellfun(@mean,times_cell_feb);

[step_freqs_mar,~] = histcounts(sample_table_mar.step_no,edges);
step_freqs_mar = step_freqs_mar';

times_cell_mar = mat2cell(sample_table_mar.t_gen,step_freqs_mar);
mean_post_mar = cellfun(@mean,times_cell_mar);

[step_freqs_apr,~] = histcounts(sample_table_apr.step_no,edges);
step_freqs_apr = step_freqs_apr';

times_cell_apr = mat2cell(sample_table_apr.t_gen,step_freqs_apr);
mean_post_apr = cellfun(@mean,times_cell_apr);

[step_freqs_may,~] = histcounts(sample_table_may.step_no,edges);
step_freqs_may = step_freqs_may';

times_cell_may = mat2cell(sample_table_may.t_gen,step_freqs_may);
mean_post_may = cellfun(@mean,times_cell_may);

[step_freqs_jun,~] = histcounts(sample_table_jun.step_no,edges);
step_freqs_jun = step_freqs_jun';

times_cell_jun = mat2cell(sample_table_jun.t_gen,step_freqs_jun);
mean_post_jun = cellfun(@mean,times_cell_jun);

[step_freqs_jul,~] = histcounts(sample_table_jul.step_no,edges);
step_freqs_jul = step_freqs_jul';

times_cell_jul = mat2cell(sample_table_jul.t_gen,step_freqs_jul);
mean_post_jul = cellfun(@mean,times_cell_jul);

[step_freqs_aug,~] = histcounts(sample_table_aug.step_no,edges);
step_freqs_aug = step_freqs_aug';

times_cell_aug = mat2cell(sample_table_aug.t_gen,step_freqs_aug);
mean_post_aug = cellfun(@mean,times_cell_aug);

figure(); hold on;
v = violinplot([mean_post_feb,mean_post_mar,mean_post_apr,mean_post_may,mean_post_jun,mean_post_jul,mean_post_aug],{'','','','','','',''},'ShowData',false,'ViolinAlpha',1);
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot; v3 = v(3).ViolinPlot; v4 = v(4).ViolinPlot; v5 = v(5).ViolinPlot; v6 = v(6).ViolinPlot; v7 = v(7).ViolinPlot;
ylabel({'Mean generation time (days)'})
l = legend([v1,v2,v3,v4,v5,v6,v7],{'February','March','April','May','June','July','August'});
ylim([0,10])

% Save comparisons

variant_comparison.data1 = mean_post_Alpha;
variant_comparison.data2 = mean_post_Delta;

variant_sd_comparison.data1 = sd_post_Alpha;
variant_sd_comparison.data2 = sd_post_Delta;

vaccine_infector_variant_comparison.data1 = mean_post_unvacc_Alpha;
vaccine_infector_variant_comparison.data2 = mean_post_onedose_Alpha;
vaccine_infector_variant_comparison.data3 = mean_post_twodose_Alpha;
vaccine_infector_variant_comparison.data4 = mean_post_unvacc_Delta;
vaccine_infector_variant_comparison.data5 = mean_post_onedose_Delta;
vaccine_infector_variant_comparison.data6 = mean_post_twodose_Delta;

vaccine_infectee_variant_comparison.data1 = mean_post_unvacc_infectee_Alpha;
vaccine_infectee_variant_comparison.data2 = mean_post_onedose_infectee_Alpha;
vaccine_infectee_variant_comparison.data3 = mean_post_twodose_infectee_Alpha;
vaccine_infectee_variant_comparison.data4 = mean_post_unvacc_infectee_Delta;
vaccine_infectee_variant_comparison.data5 = mean_post_onedose_infectee_Delta;
vaccine_infectee_variant_comparison.data6 = mean_post_twodose_infectee_Delta;

vaccine_combo_variant_comparison.data1 = mean_post_uunvacc_Alpha;
vaccine_combo_variant_comparison.data2 = mean_post_uvacc_Alpha;
vaccine_combo_variant_comparison.data3 = mean_post_vunvacc_Alpha;
vaccine_combo_variant_comparison.data4 = mean_post_vvacc_Alpha;
vaccine_combo_variant_comparison.data5 = mean_post_uu_Delta;
vaccine_combo_variant_comparison.data6 = mean_post_uv_Delta;
vaccine_combo_variant_comparison.data7 = mean_post_vu_Delta;
vaccine_combo_variant_comparison.data8 = mean_post_vv_Delta;

age_infector_variant_comparison.data1 = mean_post_age1_Alpha;
age_infector_variant_comparison.data2 = mean_post_age2_Alpha;
age_infector_variant_comparison.data3 = mean_post_age3_Alpha;
age_infector_variant_comparison.data4 = mean_post_age4_Alpha;
age_infector_variant_comparison.data5 = mean_post_age1_Delta;
age_infector_variant_comparison.data6 = mean_post_age2_Delta;
age_infector_variant_comparison.data7 = mean_post_age3_Delta;
age_infector_variant_comparison.data8 = mean_post_age4_Delta;

age_infectee_variant_comparison.data1 = mean_post_age1_infectee_Alpha;
age_infectee_variant_comparison.data2 = mean_post_age2_infectee_Alpha;
age_infectee_variant_comparison.data3 = mean_post_age3_infectee_Alpha;
age_infectee_variant_comparison.data4 = mean_post_age4_infectee_Alpha;
age_infectee_variant_comparison.data5 = mean_post_age1_infectee_Delta;
age_infectee_variant_comparison.data6 = mean_post_age2_infectee_Delta;
age_infectee_variant_comparison.data7 = mean_post_age3_infectee_Delta;
age_infectee_variant_comparison.data8 = mean_post_age4_infectee_Delta;

date_comparison.data1 = mean_post_feb;
date_comparison.data2 = mean_post_mar;
date_comparison.data3 = mean_post_apr;
date_comparison.data4 = mean_post_may;
date_comparison.data5 = mean_post_jun;
date_comparison.data6 = mean_post_jul;
date_comparison.data7 = mean_post_aug;

save('../../Results/household_gen_comp_mech1.mat','variant_comparison','variant_sd_comparison')
save('../../Results/household_gen_comp_mech2.mat','vaccine_infector_variant_comparison','vaccine_infectee_variant_comparison','vaccine_combo_variant_comparison','age_infector_variant_comparison','age_infectee_variant_comparison','date_comparison')

rmpath('../../Data')