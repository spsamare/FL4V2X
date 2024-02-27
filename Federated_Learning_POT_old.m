% Federated learning code
tic
%%
%{
gradient_global = ones(1, 2); %dont make it zero
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
%% Learning
if (mod(t_sampling, window_federated)==0)&&...
        (t_sampling>=window_federated)&&...
        (mod(t, window_sampling)==0)
    
    starting_time = 0;%t_sampling-window_federated;
    
    %% at VUEs
    active_learning_indicator = ones(total_VUE_pairs, 1);
    max_sample_global = max( max( max( maximum_queues(:, (starting_time+1):t_sampling) - threshold_Q, 0 ) ) );
    for pair = 1:total_VUE_pairs
        local_sample_size = sum( maximum_queues(pair, ...
            (starting_time+1):t_sampling)> threshold_Q );
        gradient_VUE = zeros(size(gradient_global));
        count_VUE_samples = 0;
        localParam = zeros(size(evtParam_global));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if local_sample_size>0%0 %no point of learning on empty set
            % disp(pair);
            stepsize_local = stepsize_SVRG/local_sample_size;
            local_sample_set = (maximum_queues(pair, ...
                maximum_queues(pair,(starting_time+1):t_sampling)...
                >threshold_Q) - threshold_Q)';
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
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CHANGE
            randomized_samples = randperm(local_sample_size);
            % randomized_samples = 1:local_sample_size;
            for this_sample_index = randomized_samples
                this_sample = local_sample_set(this_sample_index);
                % disp(this_sample_index);
                % disp(this_sample*localParam(2) + localParam(1)>=0);
                % disp(this_sample*evtParam_global_this_set(2) + evtParam_global_this_set(1)>=0);
                gradient_now = ...
                    Gradient_MaxLikelihood_GPD( localParam, this_sample );
                gradient_now_global = ...
                    Gradient_MaxLikelihood_GPD( evtParam_global_this_set, this_sample );
                localParam_unconstrained = localParam - stepsize_local.*( ...
                    gradient_now - gradient_now_global + gradient_global);
                % project to feasible set
                localParam_projected = Projecting_Gradient_GPD( ...
                    localParam_unconstrained, max_sample );
                % If projection changes Param, gradient needs to be calculated.
                gradient_now = (localParam - localParam_projected)./stepsize_local ...
                    + gradient_now_global - gradient_global;
                localParam = localParam_projected;
                %                         disp([this_sample_index this_sample localParam]);
                %                         disp(gradient_now);
                gradient_VUE = gradient_VUE + gradient_now;
                count_VUE_samples = count_VUE_samples + 1;
                %                         test_samples = [test_samples; this_sample];
            end %end for: randomize samples
        end %end if: non empty sample set
        if count_VUE_samples>0
            gradient_local(pair,:) = gradient_VUE/count_VUE_samples;
        else
            gradient_local(pair,:) = gradient_VUE;
            active_learning_indicator(pair) = 0;
        end
        samples_local(pair) = count_VUE_samples;
        evtParam_local(pair,:) = localParam;
    end%end for: each VUE
    
    %% at RSU
    samples_global = sum(samples_local);
    if (samples_global>0)
        gradient_global = sum( bsxfun(@times, gradient_local, samples_local), 1) ...
            / sum(samples_local);
%         gradient_global = sum( gradient_local(active_learning_indicator==1,:), 1)/...
%             samples_global;
    end
    toc;
    evtParam_global = evtParam_global + sum( ...
        repmat( samples_local(active_learning_indicator==1)/samples_global, [1, 2]).*...
        bsxfun( @plus, evtParam_local(active_learning_indicator==1,:), -evtParam_global), ...
        1)
    
    threshold_maxQ_moment_1 = threshold_Q + ...
        evtParam_global(1) / ( 1 - evtParam_global(2) );
    
    %% testing
%     Test_Plot;
%     pause(.1);
end