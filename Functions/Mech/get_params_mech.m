function params = get_params_mech(theta,params_known)

    % Recover the vector of parameters for the individual infectiousness
    % model for a given variant, params = [gamma,mu,k_inc,k_E,k_I,alpha],
    % from the vector of parameters estimated in the MCMC model fitting
    % procedure for that variant, theta = [k_E/k_inc,1/mu,alpha,beta], and
    % the vector of known parameters, params_known =
    % [k_inc,gamma,k_I,rho,x_A]
    
    k_inc = params_known(1);
    gamma = params_known(2);
    k_I = params_known(3);
    
    k_E = (theta(1))*k_inc;
    mu = 1/(theta(2));
    alpha = (theta(3));
    beta0 = (theta(4));
    
    params = [gamma,mu,k_inc,k_E,k_I,alpha,beta0];

end