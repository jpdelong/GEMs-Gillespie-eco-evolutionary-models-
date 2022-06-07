function [x_dist_out, stand_times, R_data_out, x_data_out, x_var_data_out] = GEMv2_logistic_growthalt_4(j, b, d, to_slope, bslope, dslope, cv, h_2, num_replicates, y0, t_max, Rcull, cr)

% this is a GEM for exonential growth
% the evolving trait is intrinsic rate of population growth
% call this function to run it
% NEEDS THESE FUNCTIONS: 'pick_individuals'

%% standardized time steps for storing time series
% import the time span (t_max) but decide on step lengths
stand_times = 0:1:t_max;
num_time_steps = length(stand_times); % calculate the number of standardized time steps

%% preallocate space in which to log data at standardized times
% need one for each population, one for each evolving trait, and one for
% the variance in each evolving trait
R_stand = nan(num_replicates,num_time_steps); % population size
x_stand = nan(num_replicates,num_time_steps); % trait
x_var_stand = nan(num_replicates,num_time_steps); % trait variance
x_dist_out = [];

%% OPTIONAL: keep a record of events and associated traits
compile_events_across_reps = cell(num_replicates,4);
compile_traits_across_reps = cell(num_replicates,4);

%% start Gillespie algorithm
% run for loop for each replicate simulation
for i = 1:num_replicates
    i % display replicate in the Command Window
    rng('shuffle'); % change random number seed
    
    %% preallocate for whole time series
        R = zeros(1,1e7); % 1e7 is just a large number to ensure the vector is long enough
        x_mean = nan(1e7,1); % mean trait
        x_var = nan(1e7,1); % trait variance
        t = nan(1,1e7); % time steps
        list_events = nan(1,1e7); % actual events
        list_traits = nan(1,1e7); % sequence of traits used
        
    %% specify initial conditions
        t(1) = 0; % initial time
        R(1) = y0; % initial population size
    % pull initial distribution of trait
        x_dist_init = pick_individuals(b,cv*b,R(1));
        x_dist_init(:,2) = 1; % 1 for alive, 0 for dead
        x_dist_init(:,3) = 0; % tally of births = lrs
        x_dist_init(:,4) = 0; % time of birth
        x_dist_init(:,5) = R(1); % density at time of birth
        x_dist = x_dist_init; % reset trait distribution at the start of each simulation
        x_mean(1) = median(x_dist(:,1)); % initial mean trait
        x_var(1) = var(x_dist(:,1)); % initial variance in trait
%         pHat = lognfit(x_dist(:,1));
%         x_mean(1) = exp(pHat(1));
%         x_var(1) = pHat(2);

    %% Initiate core GEM algorithm
    count = 1 % start counter to index steps while inside loop
    while t(count) < t_max
        %t(count)
        if R(count) > 0 % as long as population size is > 0, pick another individual
            whos_alive = find(x_dist(:,2)==1);
            whosnext = randi(length(whos_alive),1); % randomly choose individual from the vector
            x_next = x_dist(whos_alive(whosnext),1); % pick the trait for that individual
            list_traits(count) = x_next;
        end        
        
        d_next = to_slope*x_next^2;
%         if j == 1 % 0 variance
%             d_next = to_slope*x_next^2;
%         elseif j == 2 % variance in b only
%             d_next = to_slope*b^2;
%         elseif j == 3 % variance in d only
%             d_next = to_slope*x_next^2; % d trades off with b
%             x_next = b;
%         elseif j == 4 % variance in b and d
%             d_next = to_slope*x_next^2; % d trades off with b
%         end
        
        % set up rates of each possible event
        % birth
            b_R = max(x_next - bslope*R(count),0)*R(count);
        % death
            d_R1 = (d_next + dslope*R(count))*R(count);
        % cull
            d_R2 = cr*(max((R(count)-Rcull),0));

    % sum the events to make wheel of fortune
        CS_vector = cumsum([b_R d_R1 d_R2]);
        Slice_widths = CS_vector./CS_vector(end);
        LI = rand < Slice_widths;
        Event_index = find(LI,1,'first');
        if Event_index > 0
            list_events(count) = Event_index;
        end

    % now choose actual events
        if Event_index == 1 % choose birth of prey
            if h_2 == 0
                x_parent = (1-h_2)*mean(x_dist_init(:,1)) + h_2*x_next; % offspring trait distribution mean
            elseif h_2 > 0
                x_parent = (1-h_2)*mean(x_dist(whos_alive,1)) + h_2*x_next; % offspring trait distribution mean
            end
            off_std = sqrt(1-h_2^2)*((1-h_2)*std(x_dist_init(:,1))+h_2*sqrt(x_var(count))); % offspring trait distribution std

            offspring_start = pick_individuals(x_parent,off_std,1); % return trait
            x_dist(size(x_dist,1)+1,:) = offspring_start;%*(1+0.06*randn(1)); % return trait
            x_dist(end,2) = 1; % specify alive
            x_dist(end,3) = 0; % specify hasn't reproduced yet
            x_dist(end,4) = t(count); % specify time of birth
            x_dist(end,5) = R(count); % specify time of birth
            x_dist(whos_alive(whosnext),3) = x_dist(whos_alive(whosnext),3)+1; % increase births for this individual           
        elseif Event_index == 2 % choose death
            x_dist(whos_alive(whosnext),2) = 0; % reduce dist by lost individual
        elseif Event_index == 3 % choose cull
            x_dist(whos_alive(whosnext),2) = 0;            
        end
                
        if isempty(x_dist) == 1
            R(count+1) = 0; % need to zero out pop if x_dist is empty (otherwise = 1 for empty row vector) 
        else
            R(count+1) = sum(x_dist(:,2)); % count population size if vector is not empty
        end
        
        whos_alive = find(x_dist(:,2)==1);
        x_mean(count+1) = median(x_dist(whos_alive,1)); % calculate new mean trait
        x_var(count+1) = var(x_dist(whos_alive,1)); % calculate new variance trait
      
%         pHat = lognfit(x_dist(whos_alive,1));
%         x_mean(count+1) = exp(pHat(1));
%         x_var(count+1) = pHat(2);
        
        t(count+1) = t(count) + exp(-1/CS_vector(end))/CS_vector(end); % updating time
        count = count+1; % advance time steps
    end
        
    %% find standardized times and corresponding densities (need for ci's)    
        tmp = abs(t(~isnan(t))-stand_times');
        minima = min(tmp,[],2);
        idx = tmp==minima;
        touse = sum(idx,1);

        R_stand(i,1:length(R(touse==1))) = R(touse==1); % resource abundance at standard time
        x_stand(i,1:length(R(touse==1))) = x_mean(touse==1); % trait values at standard time
        x_var_stand(i,1:length(R(touse==1))) = x_var(touse==1); % trait values at standard time
                
        x_dist_out = [x_dist_out; x_dist];
end

% calculate ci's for time series
    upper_ci_level = 75; % choose ci levels
    lower_ci_level = 25; % choose ci levels
    R_data_out = medians_and_cis(upper_ci_level,lower_ci_level,R_stand(:,:)); % abundance      
    x_data_out = medians_and_cis(upper_ci_level,lower_ci_level,x_stand(:,:)); % trait
    x_var_data_out = medians_and_cis(upper_ci_level,lower_ci_level,x_var_stand(:,:)); % variance in trait 
