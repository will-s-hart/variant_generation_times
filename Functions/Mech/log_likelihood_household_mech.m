function l_household = log_likelihood_household_mech(f_inc,params_indiv_dir,b_cond_form,B_cond_form,mean_transmissions_form,data_struct_augmented)
        
    % Calculate the contribution from each household to the
    % log-likelihood of the augmented data, data_struct_augmented, for the
    % mechanistic model.
    
    % Throughout, arrays with the suffix "_dir" are ordered according to
    % (i) the household number assigned to each household, and (ii) the
    % (unknown) order in which household members became infected.

    % Individual parameters

    eta_indiv_dir = params_indiv_dir(:,8);
    beta_prim0_indiv_dir = params_indiv_dir(:,9);
    
    % Augmented infection and symptom onset times
    
    t_i_dir = data_struct_augmented.t_i_dir;
    t_s_dir = data_struct_augmented.t_s_dir;
    
    % Number of hosts
    
    no_hosts = length(t_i_dir);
    
    % Logical vectors indicating infection/symptom/primary status
    
    infected_dir = data_struct_augmented.infected_dir;
    uninfected_dir = ~infected_dir;
    primary_dir = data_struct_augmented.primary_dir;
    
    % Indicator matrix showing which household each individual belongs to
    
    household_indicator_mat = data_struct_augmented.household_indicator_mat;
    
    % Determine whether coprimary transmission is included
    
    coprimaries = any(beta_prim0_indiv_dir);
    
    if coprimaries
        t_start_household = round(t_i_dir(primary_dir))-0.5;
        t_start_indiv_dir = household_indicator_mat*t_start_household;
    end
    
    % Information about possible infectors
    
    poss_infectors_dir = data_struct_augmented.poss_infectors_dir;
    v = poss_infectors_dir.all; %list of all possible infectors for each host in order
    M_from = poss_infectors_dir.from_indicator_mat;
    M_to = poss_infectors_dir.to_indicator_mat;
    to_uninfected_indicator_v = ~poss_infectors_dir.to_infected_indicator;
    to_primary_indicator_v = poss_infectors_dir.to_primary_indicator;
    to_recipient_indicator_v = poss_infectors_dir.to_recipient_indicator;
    
    % Parameter values corresponding to possible infectors in v
    
    params_v = M_from*params_indiv_dir;
    params_v(to_primary_indicator_v,:) = NaN; %Not really necessary

    % Contribution to likelihood from incubation periods
        
    t_i_dir_inf = t_i_dir(infected_dir);
    t_s_dir_inf = t_s_dir(infected_dir);
        
    t_inc_inf = t_s_dir_inf-t_i_dir_inf;
    l1_inf = log(f_inc(t_inc_inf));
    
    l1 = zeros(no_hosts,1);
    l1(infected_dir) = l1_inf;
        
    % Contribution to likelihood from transmissions occurring
         
    t_tost_contribs = M_to*t_i_dir-M_from*t_s_dir; %vector of every possible generation time for each infected host
    t_inc_contribs = M_from*(t_s_dir-t_i_dir);
    
    t_tost_recipient_contribs = t_tost_contribs(to_recipient_indicator_v);
    t_inc_recipient_contribs = t_inc_contribs(to_recipient_indicator_v);
    params_v_recipient_contribs = params_v(to_recipient_indicator_v,:);
        
    t_inc_uninfected_contribs = t_inc_contribs(to_uninfected_indicator_v);
    params_v_uninfected_contribs = params_v(to_uninfected_indicator_v,:);
    
    L2a_contribs = zeros(length(v),1);
    L2a_contribs(to_recipient_indicator_v) = b_cond_form(t_tost_recipient_contribs,t_inc_recipient_contribs,params_v_recipient_contribs);
    
    L2a = eta_indiv_dir.*(M_to'*L2a_contribs);
    
    if coprimaries
        L2b = eta_indiv_dir.*beta_prim0_indiv_dir.*(t_i_dir < t_start_indiv_dir + 1);
        L2b(uninfected_dir) = 0; %not strictly necessary
    else
        L2b = zeros(no_hosts,1);
        L2b(primary_dir) = eta_indiv_dir(primary_dir)./(household_indicator_mat'*eta_indiv_dir);
    end
    
    L2 = L2a + L2b;
    L2(uninfected_dir) = 1;
    l2 = log(L2);
    
    % Contribution to likelihood from evasion of infection up to time of
    % infection (or for all time for individuals who remained uninfected)

    l3a_contribs = zeros(length(v),1);
    l3a_contribs(to_recipient_indicator_v) = B_cond_form(t_tost_recipient_contribs,t_inc_recipient_contribs,params_v_recipient_contribs);
    l3a_contribs(to_uninfected_indicator_v) = mean_transmissions_form(t_inc_uninfected_contribs,params_v_uninfected_contribs);
    
    l3a = -eta_indiv_dir.*(M_to'*l3a_contribs);
    
    if coprimaries
        l3b = -0.5*eta_indiv_dir.*beta_prim0_indiv_dir.*(1+abs(t_i_dir-t_start_indiv_dir)-abs(t_i_dir-t_start_indiv_dir-1));
        l3b(uninfected_dir) = -eta_indiv_dir(uninfected_dir).*beta_prim0_indiv_dir(uninfected_dir);
    else
        l3b = zeros(no_hosts,1);
    end
    
    l3 = l3a + l3b;
    
    % Calculate overall likelihood contribution from each individual
    
    l_indiv = l1 + l2 + l3;
    
    % Likelihood contribution from each household
    
    l_household = household_indicator_mat'*l_indiv;

    % Correction for conditioning on at least one primary
    
    if coprimaries
        l_household = l_household - log(1-exp(-household_indicator_mat'*(eta_indiv_dir.*beta_prim0_indiv_dir)));
    end
end