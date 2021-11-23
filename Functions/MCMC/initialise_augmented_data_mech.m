function data_struct_augmented = initialise_augmented_data_mech(data_struct_observed)

    % Initialise the structure array of augmented data,
    % data_struct_augmented, from the structure array of observed data,
    % data_struct_observed, for the mechanistic model
    
    % Throughout, arrays with the suffix "_dir" are ordered according to
    % (i) the household number assigned to each household, and (ii) the
    % (unknown) order in which household members became infected.
    % Otherwise, arrays are ordered according to (i) household number, and
    % (ii) the (known) order in which household members developed symptoms.
        
    % Household number that each individual belongs to, and the total
    % number of individuals
    
    household_no = data_struct_observed.household_no;
    no_hosts = length(household_no);

    household_sizes_incl = data_struct_observed.household_sizes_incl;
    household_indicator_mat = data_struct_observed.household_indicator_mat;

    % Left and right bounds for infection and symptom onset times    
    
    t_iL = data_struct_observed.t_iL;
    t_iR = data_struct_observed.t_iR;
    t_sL = data_struct_observed.t_sL;
    t_sR = data_struct_observed.t_sR;

    % Logical arrays indicating symptom status    
    
    symp = data_struct_observed.symp;
    asymp = data_struct_observed.asymp;
    
    symp_in_household = data_struct_observed.symp_in_household;
    no_symp_in_household = data_struct_observed.no_symp_in_household;
    asymp_in_household = data_struct_observed.asymp_in_household;

    % Initialise vectors t_i and t_s containing infection and symptom onset
    % times for each individual
    
    t_s = inf*ones(no_hosts,1);
    t_i = inf*ones(no_hosts,1);

    % Infection/onset times for symptomatic infected hosts
    
    t_s(symp) = t_sL(symp) + (t_sR(symp)-t_sL(symp)).*rand(sum(symp),1);
    t_i(symp) = t_s(symp) - 6;
    
    % Infection times for asymptomatic infected hosts
    
    asymp1 = asymp&(household_indicator_mat*symp_in_household);
    asymp1_in_household = (symp_in_household&asymp_in_household);
    no_asymp1_in_household = household_indicator_mat'*asymp1;
    
    r1 = rand(sum(asymp1_in_household),1);
    household_start_indices = (cumsum(household_sizes_incl)-household_sizes_incl);
    reference_hosts = household_start_indices(asymp1_in_household)+ceil(no_symp_in_household(asymp1_in_household).*r1);
    
    %MAYBE TRY TO FIND A NICER EQUIVALENT OF THIS LINE LATER
    reference_host_indiv = cell2mat(arrayfun(@(val,reps)repmat(val,reps,1),reference_hosts,no_asymp1_in_household(asymp1_in_household),'UniformOutput',false));
    
    t_i(asymp1) = t_i(reference_host_indiv)+10*(2*rand(sum(asymp1),1)-1);
    
    asymp2 = asymp&~(household_indicator_mat*symp_in_household);
    t_i(asymp2) = min(t_iR(asymp2)-10,max(t_iR(~isinf(t_iR))-10))+rand(sum(asymp2),1);
    
    % Ensure initial values lie within bounds
    
    t_i = max(min(t_i,t_iR),t_iL);
    
    % Assign values to "onset" times of asymptomatic hosts (i.e., time of
    % entry into the I stage)
    
    t_s(asymp) = t_i(asymp)+6;
    
    % Augmented data ordered by (estimated) time of infection within the
    % household
    
    [~,t_dir_host_inds] = sortrows([household_no,t_i]);
    t_i_dir = t_i(t_dir_host_inds);
    t_s_dir = t_s(t_dir_host_inds);
    
    symp_dir = symp(t_dir_host_inds);
    asymp_dir = asymp(t_dir_host_inds);

    % Store augmented data in a structure array

    data_struct_augmented = data_struct_observed;
    
    data_struct_augmented.t_i = t_i;
    data_struct_augmented.t_s = t_s;
    data_struct_augmented.t_i_dir = t_i_dir;
    data_struct_augmented.t_s_dir = t_s_dir;
    
    data_struct_augmented.t_dir_host_inds = t_dir_host_inds;
    
    data_struct_augmented.symp_dir = symp_dir;
    data_struct_augmented.asymp_dir = asymp_dir;
end