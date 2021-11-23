function mean_transmissions = mean_transmissions_form_mech(t_inc,params)
    
    % Expected infectiousness of a host, conditional on incubation period
    % t_inc, integrated over the entire course of infection, under the
    % mechanistic approach with parameters given by params (works
    % when params is either a row vector or a matrix with the same number
    % of rows as t_inc).

    params = repmat(params,size(params,1)/length(t_inc),1);
    
    gamma = params(:,1); mu = params(:,2);
    k_inc = params(:,3); k_E = params(:,4); k_I = params(:,5);
    alpha = params(:,6); beta = params(:,7);
    k_P = k_inc-k_E;
        
    mean_transmissions = gamma.*beta.*(alpha.*k_P.*mu.*t_inc+k_inc)./(alpha.*k_P.*mu+k_inc.*gamma);
end