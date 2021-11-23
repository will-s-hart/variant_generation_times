function [theta_mat,ll_vec,output] = fit_params(no_steps,steps_keep,update_theta,update_infection,update_onset,update_asymp,sample_gen_form,theta_init,data_struct_augmented_init,ll_household_init,plotting)
    
    % Use data augmentation MCMC to fit the parameters, theta, of the model
    % of infectiousness under consideration, to the household transmission
    % data.
    
    no_params_fitted = length(theta_init);
    
    % Initialise the vector of fitted parameters, theta, and the structure
    % array containing the augmented data, data_struct_augmented
    
    theta = theta_init;
    data_struct_augmented = data_struct_augmented_init;
    
    % Calculate initial likelihood contributions from each household
    
    ll_household = ll_household_init;
    
    % Vectors to hold output of fitting procedure (at steps at
    % which this is kept after burn-in and thinning)
    
    no_steps_kept = length(steps_keep);
    theta_mat = zeros(no_steps_kept,no_params_fitted); %parameters at each step
    ll_vec = zeros(no_steps_kept,1); %likelihood at each step
    
    gen_sample_table = table('Size',[no_steps_kept,5],'VariableTypes',{'cell','cell','cell','cell','cell'},'VariableNames', {'infector','infectee','t_gen','d_gen','presymp'});
    
    % Vectors to hold output of fitting procedure at every step
    
    theta_mat_all = zeros(no_steps,no_params_fitted);
    ll_vec_all = zeros(no_steps,1);
    
    % Vectors to hold information about which steps are accepted
    
    acceptance_vec = NaN*ones(no_steps,1);
    acceptance_vec_theta = NaN*ones(no_steps,1);
    acceptance_vec_infection_symp = NaN*ones(no_steps,1);
    acceptance_vec_onset = NaN*ones(no_steps,1);
    acceptance_vec_asymp = NaN*ones(no_steps,1); %updating just infection times of asymptomatic hosts in the independent transmission and symptoms model, but both infection times of asymptomatic hosts and times of entry into the I stage in the mechanistic model
    acceptance_vec_infection_asymp = NaN*ones(no_steps,1); %updating just infection times of asymptomatic hosts

    % Break steps into 100 groups to record progress
    
    step_no_mat = reshape(1:no_steps,no_steps/100,100);
    
    % If plotting enabled, set up a figure to plot output
    
    if plotting
        figure();
        subplot(2,1,1);
        subplot(2,1,2);
    end
    
    % Run chain
    
    step_no_kept = 0;
    
    for j = 1:size(step_no_mat,2)
        for i = 1:size(step_no_mat,1)

            step_no = step_no_mat(i,j);
            
            % Carry out one of 4 different steps
            
            if mod(step_no-1,4)<0.5
                
                % Update the model parameters, theta
                
                [theta,ll_household,acceptance] = update_theta(theta,data_struct_augmented,ll_household);
                
                acceptance_vec_theta(step_no) = acceptance.overall;
            elseif mod(step_no-1,4)<1.5
                
                % Update times of infection (for symptomatic infected hosts
                % in the independent transmission and symptoms model, and
                % for all infected hosts in the mechanistic model)
                
                [data_struct_augmented,ll_household,acceptance] = update_infection(theta,data_struct_augmented,ll_household);
                
                acceptance_vec_infection_symp(step_no) = acceptance.symp;
                acceptance_vec_infection_asymp(step_no) = acceptance.asymp;
            elseif mod(step_no-1,4)<2.5
                
                % Update times of symptom onset for hosts who developed
                % symptoms
                
                [data_struct_augmented,ll_household,acceptance] = update_onset(theta,data_struct_augmented,ll_household);
                
                acceptance_vec_onset(step_no) = acceptance.overall;
            else
                
                % Update augmented data for asymptomatic infected hosts
                % (just infection times for the independent transmission
                % and symptoms model, times of both infection and entry
                % into the I stage for the mechanistic model)
                
                [data_struct_augmented,ll_household,acceptance] = update_asymp(theta,data_struct_augmented,ll_household);
                
                acceptance_vec_asymp(step_no) = acceptance.overall;
                acceptance_vec_infection_asymp(step_no) = acceptance.infection;
            end
            
            % Record parameters/likelihood/acceptance from the current step
            
            theta_mat_all(step_no,:) = theta;
            ll_vec_all(step_no) = sum(ll_household);
            acceptance_vec(step_no) = acceptance.overall;
            
            % If step kept after burn-in and thinning, record further
            % output
            
            if ismember(step_no,steps_keep)
                
                step_no_kept = step_no_kept + 1;
                
                theta_mat(step_no_kept,:) = theta;
                ll_vec(step_no_kept) = sum(ll_household);
                
                % Save rng state and use different seed (ensuring different
                % seed each time) so number of steps kept does not affect
                % output
                
                s = rng;
                s1 = s; s1.Seed = s1.Seed + step_no; rng(s1);
                
                % Sample from household generation time distribution
                                
                samples = sample_gen_form(theta,data_struct_augmented);
                gen_sample_table.infector{step_no_kept} = samples.infector;
                gen_sample_table.infectee{step_no_kept} = samples.infectee;
                gen_sample_table.t_gen{step_no_kept} = samples.t_gen;
                gen_sample_table.d_gen{step_no_kept} = samples.d_gen;
                gen_sample_table.presymp{step_no_kept} = samples.presymp;
                
                % Restore rng state
                
                rng(s)
            end
        end
        
        % Display progress when an integer percentage of steps has been
        % completed
        
        fprintf('%d%% complete\n',100*step_no/no_steps);
        
        % If plotting enabled, plot output up to current step
        
        if plotting
            subplot(2,1,1); plot(1:step_no,theta_mat_all(1:step_no,:))
            subplot(2,1,2); plot(1:step_no,ll_vec_all(1:step_no))
            pause(0)
        end
    end
    
    % Calculate acceptance rates: overall, and for the different update
    % steps
    
    acceptance_vec_theta = acceptance_vec_theta(~isnan(acceptance_vec_theta));
    acceptance_vec_infection_symp = acceptance_vec_infection_symp(~isnan(acceptance_vec_infection_symp));
    acceptance_vec_asymp = acceptance_vec_asymp(~isnan(acceptance_vec_asymp));
    acceptance_vec_onset = acceptance_vec_onset(~isnan(acceptance_vec_onset));
    acceptance_vec_infection_asymp = acceptance_vec_infection_asymp(~isnan(acceptance_vec_infection_asymp));
    
    acceptance_rate_overall = mean(acceptance_vec)
    acceptance_rate_theta = mean(acceptance_vec_theta)
    acceptance_rate_infection_symp = mean(acceptance_vec_infection_symp)    
    acceptance_rate_onset = mean(acceptance_vec_onset)
    acceptance_rate_asymp = mean(acceptance_vec_asymp)
    acceptance_rate_infection_asymp = mean(acceptance_vec_infection_asymp)
        
    acceptance_rate.overall = acceptance_rate_overall;
    acceptance_rate.theta = acceptance_rate_theta;
    acceptance_rate.infection_symp = acceptance_rate_infection_symp;
    acceptance_rate.onset = acceptance_rate_onset;
    acceptance_rate.asymp = acceptance_rate_asymp;
    acceptance_rate.infection_asymp = acceptance_rate_infection_asymp;
    
    % Structure array containing additional output from the fitting
    % procedure
    
    output.acceptance_rate = acceptance_rate;
    output.data_struct_augmented_final = data_struct_augmented;
    output.gen_sample_table = gen_sample_table;
    output.theta_mat_trace = theta_mat_all(1:100:end,:);
end