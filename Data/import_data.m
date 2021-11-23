%% Import household transmission data from HOCO2 study
%% Initialisation and data import

clear all; close all; clc;

% Import data

T = readtable('Supplementary_Data.xlsx');

no_hosts_init = size(T,1);

%% Extract transmission data

% Vectors to store data

host_nos_init_incl = [];
index_host_nos_init_incl = [];
index_host_nos_new_incl = [];

household_sizes_incl = [];
household_sizes_full = [];
household_size_indiv_incl = [];
household_size_indiv_full = [];
household_no = [];
d_s = [];
d_iR = [];
infected_dir = logical([]);
symp = logical([]);
asymp = logical([]);
household_index_dates = [];

onset_diffs = [];

% Initialise values

all_processed = false;

household_no_incl = 0;

host_no_init = 0;

while ~all_processed
    
    % Loop through all households in study
    
    household_processed = false;
    
    % Vectors to hold household data
    
    household_host_nos_init = [];
    household_d_s = [];
    household_d_iR = [];
    household_symp = logical([]);
    household_asymp = logical([]);
    household_inf = logical([]);
    
    while ~household_processed
        
        % Loop through hosts in household
        
        host_no_init = host_no_init + 1;
                
        indiv_d_iR = NaN; indiv_d_s = NaN; indiv_inf = NaN; indiv_symp = NaN; indiv_asymp = NaN;
        
        % Distinguish between index cases and contacts
        
        if matches(T.Status{host_no_init},'CASE')
            
            % Index case
            
            household_index_host_no_init = host_no_init;
            
            % Definitely infected
            
            indiv_inf = true;
            indiv_included = true;
            
            % Symptom onset date (returns NaN if asymptomatic)
            
            indiv_d_s = T.SymptomOnsetDay(host_no_init);
            
            % Date of recruitment swab saved as index date for household
            % (if unavailable, take to be earliest of symptom onset date of index or study swabs)
            
            household_d_index = T.CaseRecruitmentSwabDay(host_no_init);
                        
            % Infected before symptoms (if developed) and before recruitment swab
            
            indiv_d_iR = min(indiv_d_s,household_d_index);
            
        elseif matches(T.Status{host_no_init},'CONTACT')
            
            % Contact
            
            indiv_included = true; %included for now
            
            % Symptom onset date (returns NaN if asymptomatic)
            
            indiv_d_s = T.SymptomOnsetDay(host_no_init);
            
            % Previous positive status (may actually be after HOCO swab
            % dates in a small number of cases)

            indiv_d_prev_pos = T.PreviousPositiveTestDay(host_no_init);
            indiv_prev_pos = ~isnan(indiv_d_prev_pos);
            
            if indiv_prev_pos
                
                % Days from previous positive result to household
                % recruitment
                
                prev_pos_to_index = (household_d_index-indiv_d_prev_pos);
                
                if (prev_pos_to_index > 28)
                    
                    % Don't count positive result as likely relates to
                    % previous infection
                    
                    indiv_d_prev_pos = NaN;
                    indiv_prev_pos = false;
                end
            end
            
            % Infection status and upper bound for time of infection
            % depending on test results
                        
            if matches(T.Swab1Result{host_no_init},'POSITIVE')
                
                % First swab positive so infected
                
                indiv_inf = true;
                
                % Infected before symptoms (if developed) and before 1st
                % swab
                
                indiv_d_iR = min([indiv_d_s,indiv_d_prev_pos,T.Swab1Day(host_no_init),T.Swab2Day(host_no_init),T.Swab3Day(host_no_init)]);
            
            elseif matches(T.Swab2Result{host_no_init},'POSITIVE')
                
                % Second swab positive so infected
                
                indiv_inf = true;
                
                % Infected before symptoms (if developed) and before 2nd
                % swab
                
                indiv_d_iR = min([indiv_d_s,indiv_d_prev_pos,T.Swab2Day(host_no_init),T.Swab3Day(host_no_init)]);
            
            elseif matches(T.Swab3Result{host_no_init},'POSITIVE')
                
                % Third swab positive so infected
                
                indiv_inf = true;
                                
                % Infected before symptoms (if developed) and before 3rd
                % swab
                
                indiv_d_iR = min([indiv_d_s,indiv_d_prev_pos,T.Swab3Day(host_no_init)]);
            
            elseif isempty([T.Swab1Result{host_no_init},T.Swab2Result{host_no_init},T.Swab3Result{host_no_init}])
                
                % No tests taken
                
                % Infected before symptoms, if any
                
                indiv_d_iR = min(indiv_d_s,indiv_d_prev_pos);
                
                % Infected if symptoms developed or sufficiently close
                % previous positive recorded
                
                indiv_inf = ((~isnan(indiv_d_s))|indiv_prev_pos);
                                
                % Exclude from analyses if not determined to be infected
                % via either symptoms or previous positive
                
                if ~indiv_inf
                    indiv_included = false;
                end
                
            else
                
                % At least one test taken, any/all tests negative
                
                % Infected before symptoms, if any
                
                indiv_d_iR = min(indiv_d_s,indiv_d_prev_pos);
                
                % Infected if symptoms developed or previous positive,
                % otherwise assume uninfected
                
                indiv_inf = (~isnan(indiv_d_s)|indiv_prev_pos);
                
            end
        end
        
        % Symptom status (only for infected hosts)
        
        indiv_symp = (indiv_inf & ~isnan(indiv_d_s));
        indiv_asymp = (indiv_inf & isnan(indiv_d_s));
        
        % Add data to household vectors if included
        
        if indiv_included
            
            household_host_nos_init = [household_host_nos_init;host_no_init];
            household_d_s = [household_d_s;indiv_d_s];
            household_d_iR = [household_d_iR;indiv_d_iR];
            household_symp = [household_symp;indiv_symp];
            household_asymp = [household_asymp;indiv_asymp];
            household_inf = [household_inf;indiv_inf];
        
        end
        
        % Check whether at end of household
        
        if host_no_init == no_hosts_init || matches(T.Status{host_no_init+1},'CASE')
            
            % Reached end of household
            
            household_processed = true;
            
            % Size of household and number of included household members
            
            household_size_incl = length(household_d_s);
            household_size_full = T.HouseholdSize(host_no_init);
            
            % Re-order household data, placing first symptomatic hosts
            % in order of symptom onset, then asymptomatic infected
            % hosts, then uninfected hosts

            household_data_table = table(household_d_s,household_d_iR,household_inf,household_symp,household_asymp,household_host_nos_init);
            household_data_table = sortrows(household_data_table,[1,2,3,4],{'ascend','ascend','descend','ascend'},'MissingPlacement','last');

            household_d_s = household_data_table.household_d_s;
            household_d_iR = household_data_table.household_d_iR;
            household_inf = household_data_table.household_inf;
            household_symp = household_data_table.household_symp;
            household_asymp = household_data_table.household_asymp;
            household_host_nos_init = household_data_table.household_host_nos_init;
            
            % Differences in onset dates
            
            household_onset_diffs = (diff(household_d_s(household_symp)));
            
            % Only include household if at least two included members
            
            household_included = (household_size_incl >= 2);
            
            % If household included, add data to overall data vectors
            
            if household_included
                
                household_no_incl = household_no_incl + 1;
                
                % Full and included household size
                
                household_sizes_incl = [household_sizes_incl;household_size_incl];
                household_sizes_full = [household_sizes_full;household_size_full];
                household_size_indiv_incl = [household_size_indiv_incl;repmat(household_size_incl,household_size_incl,1)];
                household_size_indiv_full = [household_size_indiv_full;repmat(household_size_full,household_size_incl,1)];
                
                % Household number for each individual
                
                household_no = [household_no;repmat(household_no_incl,household_size_incl,1)];
                
                % Numbers of household members in original data and the
                % number of the index host in both original and processed
                % data
                
                host_nos_init_incl = [host_nos_init_incl;household_host_nos_init];
                index_host_nos_init_incl = [index_host_nos_init_incl;household_index_host_no_init];
                index_host_nos_new_incl = [index_host_nos_new_incl;find(host_nos_init_incl==household_index_host_no_init)];

                % Infection and symptom data
                
                d_s = [d_s;household_d_s];
                d_iR = [d_iR;household_d_iR];
                infected_dir = [infected_dir;household_inf];
                symp = [symp;household_symp];
                asymp = [asymp;household_asymp];
                
                % Differences in symptom onset times
                
                onset_diffs = [onset_diffs;household_onset_diffs];
                
                % Recruitment date of household
                
                household_index_dates = [household_index_dates;household_d_index];
            end
        end
    end
    
    if host_no_init == no_hosts_init
        all_processed = true;
    end
end

uninfected_dir = ~infected_dir;

%% Extract additional data

no_hosts = length(d_s);
no_households = length(household_sizes_full);

% Recruitment month

recruitment_month1 = T.HouseholdRecruitmentMonth(host_nos_init_incl);

recruitment_month = zeros(no_hosts,1);
recruitment_month(matches(recruitment_month1,'January')) = 1;
recruitment_month(matches(recruitment_month1,'February')) = 2;
recruitment_month(matches(recruitment_month1,'March')) = 3;
recruitment_month(matches(recruitment_month1,'April')) = 4;
recruitment_month(matches(recruitment_month1,'May')) = 5;
recruitment_month(matches(recruitment_month1,'June')) = 6;
recruitment_month(matches(recruitment_month1,'July')) = 7;
recruitment_month(matches(recruitment_month1,'August')) = 8;
recruitment_month(matches(recruitment_month1,'September')) = 9;
recruitment_month(matches(recruitment_month1,'October')) = 10;
recruitment_month(matches(recruitment_month1,'November')) = 11;
recruitment_month(matches(recruitment_month1,'December')) = 12;

household_months = NaN(no_households,1);

for household = 1:no_households
    household_inds = find(household_no==household);        
    household_months(household) = recruitment_month(household_inds(1));
end

% Age

age_group1 = T.AgeGroup(host_nos_init_incl);

age_group = zeros(no_hosts,1);
age_group(matches(age_group1,'0-10')) = 1;
age_group(matches(age_group1,'11-18')) = 2;
age_group(matches(age_group1,'19-54')) = 3;
age_group(matches(age_group1,'55+')) = 4;

% Variant

variant1 = T.HouseholdVariant(host_nos_init_incl); %string with variant name

Alpha = matches(variant1,'Alpha');
Delta = matches(variant1,'Delta');

variant = zeros(no_hosts,1); % Set to NaN for unknown, 1 for Alpha, 2 for Delta
variant(Alpha) = 1;
variant(Delta) = 2;

household_variants = NaN(no_households,1);

for household = 1:no_households
    household_inds = find(household_no==household);        
    household_variants(household) = variant(household_inds(1));
end

% Vaccine doses and dates (note in cases where vaccine dates not available,
% sufficient time assumed to have passed for vaccine to confer full
% protection)

vaccine1 = T.VaccinationStatus(host_nos_init_incl); %string describing number of vaccine doses

unvaccinated = matches(vaccine1,'Unvaccinated');
one_dose = matches(vaccine1,'One dose');
two_dose = matches(vaccine1,'Two doses');

vaccine = NaN(no_hosts,1); % Number of doses received
vaccine(unvaccinated) = 0;
vaccine(one_dose) = 1;
vaccine(two_dose) = 2;

% Vaccine type

vaccine_type1 = T.VaccineType(host_nos_init_incl); %string describing vaccine type

vaccine_type = NaN(no_hosts,1); %1 for AZ, 2 for Pf, NaN otherwise
vaccine_type(matches(vaccine_type1,'AstraZenica')) = 1;
vaccine_type(matches(vaccine_type1,'Pfizer')) = 2;

%% Compile data into structure array

data_struct_init.household_sizes_incl = household_sizes_incl;
data_struct_init.household_sizes_full = household_sizes_full;
data_struct_init.household_size_indiv_incl = household_size_indiv_incl;
data_struct_init.household_size_indiv_full = household_size_indiv_full;
data_struct_init.household_no = household_no;
data_struct_init.d_s = d_s;
data_struct_init.d_iR = d_iR;
data_struct_init.infected_dir = infected_dir;
data_struct_init.uninfected_dir = uninfected_dir;
data_struct_init.symp = symp;
data_struct_init.asymp = asymp;

data_struct_init.household_months = household_months;
data_struct_init.recruitment_month = recruitment_month;
data_struct_init.age_group = age_group;
data_struct_init.variant = variant;
data_struct_init.vaccine = vaccine;
data_struct_init.vaccine_type = vaccine_type;

%% Save data

save('data_initial.mat','data_struct_init')