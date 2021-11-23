function p = prior_fun_mech(theta)
    
    % Calculate the prior density for the mechanistic model for a given
    % variant with estimated parameters theta

    p_p_E = betapdf(theta(1),2.1,2.1); %mean 0.5, 95% CI [0.1,0.9]
    p_mu_inv = gampdf(theta(2),7,0.7); %mean 5 days, 95% CI [2.0,9.1]
    p_alpha = gampdf(theta(3),2.65,0.75); %mean 2, 95% CI [0.4,5]
    p_beta = gampdf(theta(4),2.65,0.75); %mean 2, 95% CI [0.4,5]
    
    p = p_p_E*p_mu_inv*p_alpha*p_beta;
end