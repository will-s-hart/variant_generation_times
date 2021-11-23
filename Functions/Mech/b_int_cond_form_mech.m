function B_cond = b_int_cond_form_mech(x,t_inc,params)
    
    % Expected infectiousness of a host, conditional on incubation period
    % t_inc, integrated from times (-infinity) to x since symptom onset,
    % under the mechanistic approach with parameters given by params (works
    % when params is either a row vector or a matrix with the same number
    % of rows as x and t_inc).

    params = repmat(params,size(params,1)/length(x),1);
    
    gamma = params(:,1); mu = params(:,2);
    k_inc = params(:,3); k_E = params(:,4); k_I = params(:,5);
    alpha = params(:,6); beta = params(:,7);
    k_P = k_inc-k_E;
    
    C = k_inc.*gamma.*mu./(alpha.*k_P.*mu+k_inc.*gamma);
        
    ind_m = (x<0);
    ind_p = ~ind_m;
    
    x_m = x(ind_m);
    t_inc_m = t_inc(ind_m);
    beta_m = beta(ind_m);
    alpha_m = alpha(ind_m);
    C_m = C(ind_m);
    k_P_m = k_P(ind_m);
    k_E_m = k_E(ind_m);
    k_inc_m = k_inc(ind_m);
    
    x_p = x(ind_p);
    t_inc_p = t_inc(ind_p);
    beta_p = beta(ind_p);
    alpha_p = alpha(ind_p);
    C_p = C(ind_p);
    k_P_p = k_P(ind_p);
    k_inc_p = k_inc(ind_p);
    k_I_p = k_I(ind_p);
    mu_p = mu(ind_p);
    
    f_m1 = alpha_m.*C_m.*beta_m.*x_m.*(1-betacdf(-x_m./t_inc_m,k_P_m,k_E_m));
    f_m2 = (alpha_m.*C_m.*beta_m.*k_P_m.*t_inc_m./k_inc_m).*(1-betacdf(-x_m./t_inc_m,k_P_m+1,k_E_m));
    f_m = f_m1+f_m2;
    
    f_p1 = C_p.*beta_p.*alpha_p.*k_P_p.*t_inc_p./k_inc_p;
    f_p2 = C_p.*beta_p.*x_p.*(1-gamcdf(x_p,k_I_p,1./(k_I_p.*mu_p)));
    f_p3 = (C_p.*beta_p./mu_p).*gamcdf(x_p,k_I_p+1,1./(k_I_p.*mu_p));
    f_p = f_p1+f_p2+f_p3;
    
    B_cond = zeros(size(x));
    B_cond(ind_m) = f_m;
    B_cond(ind_p) = f_p;
end