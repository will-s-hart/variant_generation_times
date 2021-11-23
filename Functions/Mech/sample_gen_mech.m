function out = sample_gen_mech(params_indiv_dir,b_cond_form,data_struct_augmented)
    
    % Sample realised generation times for the mechanistic model with
    % augmented data data_struct_augmented.
    
    % Throughout, arrays with the suffix "_dir" are ordered according to
    % (i) the household number assigned to each household, and (ii) the
    % (unknown) order in which household members became infected.
    
    % Individual parameters
    
    beta_prim0_indiv_dir = params_indiv_dir(:,9);
    
    % Augmented infection and symptom onset times
    
    t_i_dir = data_struct_augmented.t_i_dir;
    d_i_dir = round(t_i_dir);
    t_s_dir = data_struct_augmented.t_s_dir;
        
    % Number of hosts
    
    no_hosts = length(t_i_dir);
    
    % Logical vectors indicating symptom status
    
    primary_dir = data_struct_augmented.primary_dir;
    
    % Indicator matrix showing which household each individual belongs to
    
    household_indicator_mat = data_struct_augmented.household_indicator_mat;
    
    % Information about possible infectors
    
    poss_infectors_dir = data_struct_augmented.poss_infectors_dir;
    v = poss_infectors_dir.all; %list of all possible infectors for each host in order
    v_to = poss_infectors_dir.to;
    M_from = poss_infectors_dir.from_indicator_mat;
    M_to = poss_infectors_dir.to_indicator_mat;
    to_primary_indicator_v = poss_infectors_dir.to_primary_indicator;
    to_recipient_indicator_v = poss_infectors_dir.to_recipient_indicator;
    
    % Parameter values corresponding to possible infectors in v
    
    params_v = M_from*params_indiv_dir;
    params_v(to_primary_indicator_v) = NaN; %Not really necessary
    
    % Vectors of every possible generation time for each infectee (with
    % entries corresponding to infection by every possible infector), and
    % of the incubation period (which is infinite for asymptomatic hosts)
    % of the corresponding infectee.
    
    t_gen_contribs = (M_to-M_from)*t_i_dir; %vector of every possible generation time for each infected host
    d_gen_contribs = (M_to-M_from)*d_i_dir;
    t_inc_contribs = M_from*(t_s_dir-t_i_dir);
    
    % Calculate likelihood contributions corresponding to infection of each
    % individual    
    
    t_gen_recipient_contribs = t_gen_contribs(to_recipient_indicator_v);
    d_gen_recipient_contribs = d_gen_contribs(to_recipient_indicator_v);
    t_inc_recipient_contribs = t_inc_contribs(to_recipient_indicator_v);
    params_v_recipient_contribs = params_v(to_recipient_indicator_v,:);

    L2ai_contribs = zeros(length(v),1);
    L2ai_contribs(to_recipient_indicator_v) = b_cond_form(t_gen_recipient_contribs-t_inc_recipient_contribs,t_inc_recipient_contribs,params_v_recipient_contribs);
    L2ai = M_to'*L2ai_contribs;
    
    if any(beta_prim0_indiv_dir)
    
        % Primary correction

        t_start_household = round(t_i_dir(primary_dir))-0.5;
        t_start_indiv_dir = household_indicator_mat*t_start_household;

        L2aii = beta_prim0_indiv_dir.*(t_i_dir < t_start_indiv_dir + 1);    

        % Sample hosts who are infected from within the household (i.e., not
        % co-primaries)

        L2a = L2ai + L2aii;

        recipient1_dir =  (rand(no_hosts,1).*L2a < L2ai);
        recipients1_dir = find(recipient1_dir);

        no_recipients1 = length(recipients1_dir);

        to_recipient1_indicator_v = arrayfun(@(x)any(ismembertol(recipients1_dir,x)),v_to);

        t_gen_recipient1_contribs = t_gen_contribs(to_recipient1_indicator_v);
        d_gen_recipient1_contribs = d_gen_contribs(to_recipient1_indicator_v);
        t_inc_recipient1_contribs = t_inc_contribs(to_recipient1_indicator_v);        
    else
        no_recipients1 = sum(data_struct_augmented.recipient_dir);
        to_recipient1_indicator_v = to_recipient_indicator_v;
        t_gen_recipient1_contribs = t_gen_recipient_contribs;
        d_gen_recipient1_contribs = d_gen_recipient_contribs;
        t_inc_recipient1_contribs = t_inc_recipient_contribs;
    end
    
    if no_recipients1 > 0
                
        % Use likelihood contributions to calculate weights indicating the
        % relative probabilities of infection by different possible infectiors    

        weights_all = L2ai_contribs./(M_to*L2ai);
        weights_gen = weights_all(to_recipient1_indicator_v);

        % Sample from generation time distribution

        r = rand(no_recipients1,1)+(0:(no_recipients1-1))';
        inds = discretize(r,[0;cumsum(weights_gen)]);
        
        t_gen_sample = t_gen_recipient1_contribs(inds);
        d_gen_sample = d_gen_recipient1_contribs(inds);
                
        % Corresponding infectors and infectees, in estimated order of
        % infection
        
        v_recipient1_contribs = v(to_recipient1_indicator_v);
        v_to_recipient1_contribs = v_to(to_recipient1_indicator_v);

        infector_dir_sample = v_recipient1_contribs(inds);
        infectee_dir_sample = v_to_recipient1_contribs(inds);
        
        % Infectors and infectees in original ordering

        t_dir_host_inds = data_struct_augmented.t_dir_host_inds;

        infector_sample = t_dir_host_inds(infector_dir_sample); %hard to get head around but this correct
        infectee_sample = t_dir_host_inds(infectee_dir_sample);
        
        % Determine whether/not transmissions are presymptomatic (among
        % those from infectors developing symptoms
                
        t_inc_sample = t_inc_recipient1_contribs(inds);
        presymp_sample = double((t_gen_sample < t_inc_sample));
        
        asymp_sample = data_struct_augmented.asymp_dir(infector_dir_sample);
        presymp_sample(asymp_sample) = NaN;

    else
        
        infector_sample = [];
        infectee_sample = [];
        t_gen_sample = [];
        d_gen_sample = [];
    end 
    
    
    % Output
    
    out.infector = infector_sample;
    out.infectee = infectee_sample;
    out.t_gen = t_gen_sample;
    out.d_gen = d_gen_sample;
    out.presymp = presymp_sample;
end