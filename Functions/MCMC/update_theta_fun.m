function [theta_new,ll_household_new,acceptance] = update_theta_fun(theta_old,data_struct_augmented,ll_household_old,ll_household_form,theta_prop_cov_mat,prior_fun)

    % Update the vector of fitted model parameters, theta, in the parameter
    % fitting procedure
    
    % Propose new values of each entry of theta using a multivariate normal
    % proposal distribution
    
    no_params_fitted = length(theta_old);
    theta_prop = theta_old + mvnrnd(zeros(1,no_params_fitted),theta_prop_cov_mat);

    % Calculate the loglikelihood using the proposed parameters
    
    ll_household_prop = ll_household_form(theta_prop,data_struct_augmented);

    % Calculate the logarithm of the ratio between the proposed and
    % previous likelihoods
    
    a1 = exp(sum(ll_household_prop)-sum(ll_household_old));
    
    % Calculate the acceptance probability, a, for the new parameters,
    % accounting for priors
    
    a = a1*prior_fun(theta_prop)/prior_fun(theta_old);
    
    % Accept the proposed parameters and onset times with
    % probability a
    
    if rand < a

        % Accept the proposed parameters
        
        theta_new = theta_prop;
        
        % New likelihood
        
        ll_household_new = ll_household_prop;
        
        % Information about acceptance
        
        acceptance.overall = 1;
    else
        
        % Reject the proposed parameters
        
        theta_new = theta_old;
        
        % New likelihood
        
        ll_household_new = ll_household_old;
        
        % Information about acceptance
        
        acceptance.overall = 0;
    end
end