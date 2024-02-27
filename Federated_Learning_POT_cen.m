% Federated learning code
% tic
%%
%{
gradient_global = [1 1000]; %dont make it zero
evtParam_global = [50 0]; % [scale shape]
samples_global = 0;
%local models
gradient_local = repmat( gradient_global, [total_VUE_pairs, 1]); %
evtParam_local = repmat( evtParam_global, [total_VUE_pairs, 1]); % [location scale shape]
%}
%{
test_samples = [];%************????????????????
if (t==window_federated*window_sampling)
    these_samples = max( maximum_queues(:, 1:t_sampling), [], 1);
    % these_samples = sort(these_samples);
    pd = fitdist( these_samples', 'GeneralizedExtremeValue');
    disp(pd);
    evtParam_global = [pd.mu pd.sigma pd.k];
    evtParam_local = repmat( evtParam_global, [total_VUE_pairs, 1]);
end
%}
%% Save and set back the randomize state
current_random_state = rng;
rng(t_sampling);

%% Learning
if (mod(t_sampling, window_federated)==0)&&...
        (t_sampling>=window_federated)&&...
        (mod(t, window_sampling)==0)
    
    starting_time = 0;%
    %     starting_time = t_sampling-window_federated;
    
    %% at VUEs
    active_learning_indicator = 1;
    max_sample_global = max( max_sample_global, ...
        max( max( max( maximum_queues(:, ...
        (starting_time+1):t_sampling) - threshold_Q, 0 ) ) ) + 1);
    %     for pair = 1:total_VUE_pairs
    all_samples =  max( maximum_queues(:, ...
        (starting_time+1):t_sampling) - threshold_Q, 0);
    all_samples = all_samples(all_samples>0);
    local_sample_size = length(all_samples);
    gradient_VUE = zeros(size(gradient_global));
    count_VUE_samples = 0;
    localParam = zeros(size(evtParam_global));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if local_sample_size>0%0 %no point of learning on empty set
        % disp(pair);
        stepsize_local = stepsize_SVRG/local_sample_size;
        local_sample_set = all_samples; %%%%%%?????
        max_sample = max_sample_global;
        %             max_sample = max(local_sample_set);
        %             max_sample = local_sample_set;
        % Adjust parameters to fit into local data set
        % otherwise gradient becomes complex number
        evtParam_global_this_set = Projecting_Gradient_GPD( ...
            evtParam_global, max_sample );
        % localParam = evtParam_global;
        localParam = evtParam_global_this_set - stepsize_local.* ...
            gradient_global;
        localParam = Projecting_Gradient_GPD( ...
            localParam, max_sample );
        % SVGR needs couple of iterations over the same sample set for accurate
        % output.
        gradient_global_this = gradient_global; %loca version of global one
        for this_iteration = 1:iterations_SVGR
            gradient_this = zeros(size(gradient_global));
            count_VUE_samples = 0;
            randomized_samples = randperm(local_sample_size);
            for this_sample_index = randomized_samples
                %-------------
                counter_federated = counter_federated + 1;
                %-------------
                this_sample = local_sample_set(this_sample_index);
                gradient_now = ...
                    Gradient_MaxLikelihood_GPD( localParam, this_sample );
                gradient_now_global = ...
                    Gradient_MaxLikelihood_GPD( evtParam_global_this_set, this_sample );
                localParam_unconstrained = localParam - stepsize_local.*( ...
                    gradient_now - gradient_now_global + gradient_global_this);
                % project to feasible set
                localParam_projected = Projecting_Gradient_GPD( ...
                    localParam_unconstrained, max_sample );
                %                    % If projection changes Param, gradient needs to be calculated.
                %                     gradient_now = (localParam - localParam_projected)./stepsize_local ...
                %                         + gradient_now_global - gradient_global_this;
                localParam = localParam_projected;
                %                         disp([this_sample_index this_sample localParam]);
                %                         disp(gradient_now);
                gradient_this = gradient_this + gradient_now;
                count_VUE_samples = count_VUE_samples + 1;
                %                         test_samples = [test_samples; this_sample];
            end %end for: randomize samples
            gradient_global_this = gradient_this / count_VUE_samples;
            gradient_VUE = (1 - 1/this_iteration) * gradient_VUE + ...
                gradient_global_this / this_iteration;
            evtParam_global_this_set = localParam;
            
        end %end for: iterations for SVGR
        
    end %end if: non empty sample set
    if count_VUE_samples>0
        gradient_local = gradient_VUE/count_VUE_samples;
    else
        gradient_local= gradient_VUE;
        active_learning_indicator = 0;
    end
    
    samples_global = count_VUE_samples;
    gradient_global = gradient_local;
    %     toc;
    evtParam_global = localParam;
    
    threshold_maxQ_moment_1 = threshold_Q + ...
        evtParam_global(1) / ( 1 - evtParam_global(2) );
    
    %% testing
    %{
    these_samples = max( maximum_queues(:, ...
        (starting_time+1):t_sampling) - threshold_Q, 0 );
    these_samples = these_samples( these_samples > 0 );
    disp( Test_MLE( evtParam_global(1), evtParam_global(2)) );
    %     Test_Plot;
    %     pause(.1);
    %}
end

%% restore original random state
rng(current_random_state);