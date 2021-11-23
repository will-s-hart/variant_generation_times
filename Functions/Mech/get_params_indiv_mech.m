function params_indiv_dir = get_params_indiv_mech(theta,params_known,data_struct_augmented)
    
    % Calculate a matrix whose rows give the individual parameters
    % corresponding to each individual (where individuals are ordered first
    % by household, and within each household by imputed infection time).
    % The entries of each row give the following parameters for the
    % individual: [gamma,mu,k_inc,k_E,k_I,alpha,beta,eta,beta_prim0]. Here,
    % we do not account for co-primary cases, so beta_prim0 = 0 for all
    % individuals. Out of the remaining model parameters, all depend on
    % variant, while the overall transmissibility beta = beta0/(n-1) also
    % depends on household size (n) and whether the individual develops
    % symptoms, and the relative susceptibility eta depends on vaccination
    % status.
    
    % Data
    
    household_size_indiv_full_dir = data_struct_augmented.household_size_indiv_full;
    uninfected_dir = ~data_struct_augmented.infected_dir;
    asymp_dir = data_struct_augmented.asymp_dir;
    
    no_hosts = length(household_size_indiv_full_dir);
    
    variant_dir = data_struct_augmented.variant;
    
    Alpha_dir = (variant_dir==1);
    Delta_dir = (variant_dir==2);
    
    t_dir_host_inds = data_struct_augmented.t_dir_host_inds;
    
    vaccine_dir = data_struct_augmented.vaccine(t_dir_host_inds);
    unvaccinated_dir = (vaccine_dir == 0);
    one_dose_dir = (vaccine_dir == 1);
    two_dose_dir = (vaccine_dir == 2);

    vaccine_type_dir = data_struct_augmented.vaccine_type(t_dir_host_inds);
    AZ_dir = (vaccine_type_dir == 1);
    Pf_dir = (vaccine_type_dir == 2);
    
    % Assumed parameters
    
    k_inc = params_known(1);
    gamma = params_known(2);
    k_I = params_known(3);
    rho = params_known(4);
    x_A = params_known(5);
    
    beta_prim0 = 0;
    
    eta1_AZ_Alpha = 1-0.63;
    eta2_AZ_Alpha = 1-0.79;
    
    eta1_Pf_Alpha = 1-0.59;
    eta2_Pf_Alpha = 1-0.78;
    
    eta1_AZ_Delta = 1-0.46;
    eta2_AZ_Delta = 1-0.67;
    
    eta1_Pf_Delta = 1-0.57;
    eta2_Pf_Delta = 1-0.80;
    
    % Estimated parameters
    
    k_E1 = (theta(1))*k_inc;
    mu1 = 1/(theta(2));
    alpha1 = (theta(3));
    beta01 = (theta(4));
    
    k_E2 = (theta(5))*k_inc;
    mu2 = 1/(theta(6));
    alpha2 = (theta(7));
    beta02 = (theta(8));
    
    % Individual parameters
    
    gamma_indiv_dir = repmat(gamma,no_hosts,1);
    gamma_indiv_dir(uninfected_dir) = NaN;
        
    mu_indiv_dir = NaN(no_hosts,1);
    mu_indiv_dir(Alpha_dir) = mu1;
    mu_indiv_dir(Delta_dir) = mu2;
    mu_indiv_dir(uninfected_dir) = NaN;
        
    k_inc_indiv_dir = repmat(k_inc,no_hosts,1);
    k_inc_indiv_dir(uninfected_dir) = NaN;
        
    k_E_indiv_dir = NaN(no_hosts,1);
    k_E_indiv_dir(Alpha_dir) = k_E1;
    k_E_indiv_dir(Delta_dir) = k_E2;
    k_E_indiv_dir(uninfected_dir) = NaN;
        
    k_I_indiv_dir = repmat(k_I,no_hosts,1);
    k_I_indiv_dir(uninfected_dir) = NaN;
        
    alpha_indiv_dir = NaN(no_hosts,1);
    alpha_indiv_dir(Alpha_dir) = alpha1;
    alpha_indiv_dir(Delta_dir) = alpha2;
    alpha_indiv_dir(uninfected_dir) = NaN;
    
    % Overall infectiousness
    
    beta0_indiv_dir = NaN(no_hosts,1);
    beta0_indiv_dir(Alpha_dir) = beta01;
    beta0_indiv_dir(Delta_dir) = beta02;
    beta0_indiv_dir(uninfected_dir) = NaN;
    beta0_indiv_dir(asymp_dir) = x_A*beta0_indiv_dir(asymp_dir);
    beta_indiv = beta0_indiv_dir./((household_size_indiv_full_dir-1).^rho);
    
    % Relative susceptibility
    
    eta_indiv_dir = NaN(no_hosts,1);
    
    eta_indiv_dir(unvaccinated_dir) = 1;
    
    eta_indiv_dir(one_dose_dir&AZ_dir&Alpha_dir) = eta1_AZ_Alpha;
    eta_indiv_dir(two_dose_dir&AZ_dir&Alpha_dir) = eta2_AZ_Alpha;

    eta_indiv_dir(one_dose_dir&Pf_dir&Alpha_dir) = eta1_Pf_Alpha;
    eta_indiv_dir(two_dose_dir&Pf_dir&Alpha_dir) = eta2_Pf_Alpha;

    eta_indiv_dir(one_dose_dir&AZ_dir&Delta_dir) = eta1_AZ_Delta;
    eta_indiv_dir(two_dose_dir&AZ_dir&Delta_dir) = eta2_AZ_Delta;

    eta_indiv_dir(one_dose_dir&Pf_dir&Delta_dir) = eta1_Pf_Delta;
    eta_indiv_dir(two_dose_dir&Pf_dir&Delta_dir) = eta2_Pf_Delta;
    
    % Co-primary transmission rate
    
    beta_prim0_indiv_dir = repmat(beta_prim0,no_hosts,1);
    
    % Matrix of individual parameters
    
    params_indiv_dir = [gamma_indiv_dir,mu_indiv_dir,k_inc_indiv_dir,k_E_indiv_dir,k_I_indiv_dir,alpha_indiv_dir,beta_indiv,eta_indiv_dir,beta_prim0_indiv_dir];
end