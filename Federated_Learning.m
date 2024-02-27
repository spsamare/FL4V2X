% Federated learning code

%%
test_samples = [];%************????????????????
if (t==window_federated*window_sampling)
    these_samples = max( maximum_queues(:, 1:t_sampling), [], 1);
    % these_samples = sort(these_samples);
    pd = fitdist( these_samples', 'GeneralizedExtremeValue');
    disp(pd);
    evtParam_global = [pd.mu pd.sigma pd.k];
    evtParam_local = repmat( evtParam_global, [total_VUE_pairs, 1]);
end
%% Learning
if (mod(t_sampling, window_federated)==0)&&(t_sampling>window_federated)
    %% at VUEs
    active_learning_indicator = ones(total_VUE_pairs, 1);
    for pair = 1:total_VUE_pairs
        local_sample_size = sum( maximum_queues(pair, 1:t_sampling)>0 );
        gradient_VUE = zeros(size(gradient_global));
        count_VUE_samples = 0;
        localParam = zeros(size(evtParam_global));
        if local_sample_size>0 %no point of learning on empty set
            disp(pair);
            stepsize_local = stepsize_SVRG/local_sample_size;
            % localParam = evtParam_global;
            localParam = evtParam_global - stepsize_local*( ...
                + gradient_global);
            local_sample_set = (maximum_queues(pair, ...
                maximum_queues(pair,:)>0))';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CHANGE
            % randomized_samples = randperm(local_sample_size);
            randomized_samples = 1:local_sample_size;
            for this_sample_index = randomized_samples
                this_sample = local_sample_set(this_sample_index);
                temp = localParam;
                %present NaN-Inf issues
                if ( 1 + localParam(3)*(this_sample...
                        -localParam(1))/localParam(2) > 0 ) && ...
                        ( 1 + evtParam_global(3)*(this_sample...
                        -evtParam_global(1))/evtParam_global(2) > 0 )
                    gradient_now = ...
                        Gradient_MaxLikelihood_GEVD( localParam, this_sample );
                    localParam_temp = localParam - stepsize_local*( ...
                        gradient_now - ...
                        Gradient_MaxLikelihood_GEVD( evtParam_global, this_sample )...
                        + gradient_global);
                    localParam_temp = Projecting_Gradient( localParam_temp, ...
                        this_sample ); % project to feasible set
                    if ( 1 + localParam_temp(3)*(this_sample...
                            -localParam_temp(1))/localParam_temp(2) > 0 )
                        localParam = localParam_temp;
                        disp([this_sample_index this_sample localParam]);
                        disp(gradient_now);
                        gradient_VUE = gradient_VUE + gradient_now;
                        count_VUE_samples = count_VUE_samples + 1;
                        test_samples = [test_samples; this_sample];
                    else
                        disp([this_sample_index this_sample localParam_temp]);
                        disp(gradient_now);
                        disp('Rejected');
                    end
                end %end if: gradient update
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
        gradient_global = sum( gradient_local(active_learning_indicator==1,:), 1)/...
            samples_global;
    end
    evtParam_global = evtParam_global + sum( ...
        repmat( samples_local(active_learning_indicator==1)/samples_global, [1, 3]).*...
        bsxfun( @plus, evtParam_local(active_learning_indicator==1,:), -evtParam_global), ...
        1);
end