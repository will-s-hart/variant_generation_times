% Assumed values of model parameters (other than those parameter values
% that were fitted to the data)

% Lognormal incubation period from McAloon et al

inc_mu = 1.63;
inc_sigma = 0.5;
f_inc_logn = @(t_inc)lognpdf(t_inc,inc_mu,inc_sigma);
inc_mean = exp(inc_mu+0.5*inc_sigma^2);
inc_var = (exp(inc_sigma^2)-1)*exp(2*inc_mu+inc_sigma^2);

% Gamma incubation period with the same mean and standard deviation

inc_shape = inc_mean^2/inc_var;
inc_scale = inc_var/inc_mean;
f_inc_gam = @(t_inc)gampdf(t_inc,inc_shape,inc_scale);

% Assumed parameter values common to both models

x_A = 0.35; %relative infectiousness of asymptomatic infected hosts
rho = 1; %transmission scales with (household size-1)^(-rho)

% Assumed parameter values in mechanistic model

k_inc = inc_shape;
gamma = 1/(k_inc*inc_scale);
k_I = 1;

params_known = [k_inc,gamma,k_I,rho,x_A]; %vector of known parameters

% Save to .mat file
save('assumed_parameters.mat','inc_shape','inc_scale','inc_mu','inc_sigma','f_inc_logn','f_inc_gam','params_known','k_inc','gamma','k_I','rho','x_A')