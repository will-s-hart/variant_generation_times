% Format imported data into a structure array.

% Throughout, arrays with the suffix "_dir" are ordered according to (i)
% the household number assigned to each household, and (ii) the (unknown)
% order in which household members became infected. Otherwise, arrays are
% ordered according to (i) household number, and (ii) the (known) order in
% which household members developed symptoms.

clear all; close all; clc;

% Load imported data in file "data_initial.mat" (this file is created by
% running "import_data.m")

load('data_initial.mat','data_struct_init')

household_sizes_incl = data_struct_init.household_sizes_incl;
household_sizes_full = data_struct_init.household_sizes_full;
household_size_indiv_incl = data_struct_init.household_size_indiv_incl;
household_size_indiv_full = data_struct_init.household_size_indiv_full;
household_no = data_struct_init.household_no;
d_s = data_struct_init.d_s;
d_iR = data_struct_init.d_iR;
infected_dir = data_struct_init.infected_dir;
uninfected_dir = data_struct_init.uninfected_dir;
symp = data_struct_init.symp;
asymp = data_struct_init.asymp;

household_months = data_struct_init.household_months;
recruitment_month = data_struct_init.recruitment_month;
age_group = data_struct_init.age_group;
variant = data_struct_init.variant;
vaccine = data_struct_init.vaccine;
vaccine_type = data_struct_init.vaccine_type;

% Number of hosts and households

no_hosts = length(d_s);
no_households = length(household_sizes_incl);

% Set onset days/infection day bounds to infinity for hosts who did not
% develop symptoms/become infected

d_s(isnan(d_s)) = inf;
d_iR(isnan(d_iR)) = inf;

% Lower and upper bounds for symptom onset and infection times

t_sL = d_s-0.5;
t_sR = d_s+0.5;

t_iL = -inf*ones(no_hosts,1);
t_iR = d_iR+0.5;

% Indicator matrix showing which household each individual belongs to

blocks = mat2cell(sparse(ones(no_hosts,1)),household_sizes_incl,1);
household_indicator_mat = blkdiag(blocks{:});

% Calculate the number of infected and symptomatic/asymptomatic hosts in
% each household, and working in the (unknown) order of infection,
% determine the possible (household) infectors for each individual (i.e.,
% household members who were infected before that individual).

no_infected_in_household = zeros(no_households,1);
no_symp_in_household = zeros(no_households,1);
no_asymp_in_household = zeros(no_households,1);
no_poss_infectors_dir = zeros(no_hosts,1);
poss_infectors_dir_cell = cell(no_hosts,1); %cell array to populate with indices of possible infectors for each individual
primary_dir = false(no_hosts,1); %logical array to indicate which individual is the primary case in the household

for i = 1:no_households
    
    in_household = (household_no==i);
    infected_hosts_in_household = find(in_household.*infected_dir);
    symp_hosts_in_household = find(in_household.*symp);
    asymp_hosts_in_household = find(in_household.*asymp);
    uninfected_hosts_in_household = find(in_household.*uninfected_dir);
    
    no_infected_in_household(i) = length(infected_hosts_in_household);
    no_symp_in_household(i) = length(symp_hosts_in_household);
    no_asymp_in_household(i) = length(asymp_hosts_in_household);
    
    no_poss_infectors_dir(infected_hosts_in_household(1)) = 1;
    poss_infectors_dir_cell{infected_hosts_in_household(1)} = 0; %use 0 to denote infection from outside the household
    primary_dir(infected_hosts_in_household(1)) = true;
    
    for j = 2:length(infected_hosts_in_household)
        poss_infectors_dir_cell{infected_hosts_in_household(j)} = infected_hosts_in_household(1:(j-1));
        no_poss_infectors_dir(infected_hosts_in_household(j)) = length(poss_infectors_dir_cell{infected_hosts_in_household(j)});
    end
    
    for j = 1:length(uninfected_hosts_in_household)
        poss_infectors_dir_cell{uninfected_hosts_in_household(j)} = infected_hosts_in_household;
        no_poss_infectors_dir(uninfected_hosts_in_household(j)) = length(infected_hosts_in_household);
    end
end

% Infected individuals who were not primaries (working in the unknown
% order of infection)

recipient_dir = infected_dir&(~primary_dir);

% Logical arrays indicating whether or not each household contained any
% symptomatic or asymptomatic infected hosts

symp_in_household = any(no_symp_in_household,2);
asymp_in_household = any(no_asymp_in_household,2);

% Create a vector v, enumerating each values within the cell array listing
% the possible infectors for each individual (in the order of infection)

v = cell2mat(poss_infectors_dir_cell);

% Indicator matrix giving the index of the potential infector corresponding
% to each entry of v

M1 = spalloc(length(v),no_hosts,length(v)-no_households);

for i = 1:length(v)
    
    j = v(i);
    
    if j > 0
        M1(i,j) = 1;
    end
end

% Indicator matrix giving the index of the potential infectee corresponding
% to each entry of v

blocks = mat2cell(sparse(ones(length(v),1)),no_poss_infectors_dir,1);
M2 = blkdiag(blocks{:});

[~,v_to] = ind2sub(size(M2),find(M2));

% Logical vectors indicating, respectively, whether or not each entry of v
% corresponds to: (i) the potential infector of a host who actually became
% infected; (ii) the infection of the primary case (denoted by a 0 in v);
% (iii) the potential infector of a non-primary case who became infected

v_to_infected_indicator = logical(M2*infected_dir);
v_to_primary_indicator = (v==0);
v_to_recipient_indicator = (v_to_infected_indicator&~v_to_primary_indicator);

% The size of the household that each entry of v corresponds to

household_size_v = M2*household_size_indiv_full;

% Structure array containing information about possible infectors

poss_infectors_dir.cell = poss_infectors_dir_cell;
poss_infectors_dir.all = v;
poss_infectors_dir.from_indicator_mat = M1;
poss_infectors_dir.to_indicator_mat = M2;
poss_infectors_dir.to_primary_indicator = v_to_primary_indicator;
poss_infectors_dir.to_infected_indicator = v_to_infected_indicator;
poss_infectors_dir.to_recipient_indicator = v_to_recipient_indicator;
poss_infectors_dir.household_size = household_size_v;

poss_infectors_dir.to = v_to;

% Structure array containing all observed data

data_struct_observed.t_iL = t_iL;
data_struct_observed.t_iR = t_iR;
data_struct_observed.t_sL = t_sL;
data_struct_observed.t_sR = t_sR;

data_struct_observed.household_sizes_incl = household_sizes_incl;
data_struct_observed.household_sizes_full = household_sizes_full;
data_struct_observed.household_size_indiv_incl = household_size_indiv_incl;
data_struct_observed.household_size_indiv_full = household_size_indiv_full;
data_struct_observed.household_no = household_no;
data_struct_observed.household_indicator_mat = household_indicator_mat;
data_struct_observed.no_infected_in_household = no_infected_in_household;
data_struct_observed.no_symp_in_household = no_symp_in_household;
data_struct_observed.symp_in_household = symp_in_household;
data_struct_observed.no_asymp_in_household = no_asymp_in_household;
data_struct_observed.asymp_in_household = asymp_in_household;
data_struct_observed.primary_dir = primary_dir;
data_struct_observed.infected_dir = infected_dir;
data_struct_observed.recipient_dir = recipient_dir;
data_struct_observed.symp = symp;
data_struct_observed.asymp = asymp;
data_struct_observed.no_poss_infectors_dir = no_poss_infectors_dir;
data_struct_observed.poss_infectors_dir = poss_infectors_dir;

data_struct_observed.recruitment_month = recruitment_month;
data_struct_observed.household_months = household_months;
data_struct_observed.age_group = age_group;
data_struct_observed.variant = variant;
data_struct_observed.vaccine = vaccine;
data_struct_observed.vaccine_type = vaccine_type;

% Save data to .mat file

save('data.mat','data_struct_observed')