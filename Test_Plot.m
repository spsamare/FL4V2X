% Test plots
% close all

%{
%% reliability
max( mean(Queues(:,1:t)>threshold_Q,2) - violation_probability, 0 )
figure()
cdfplot( 1 - mean(Queues(:,1:7000)>threshold_Q,2))
hold on
line([1-violation_probability 1-violation_probability], [0 1])


%% Transmit powers

power_temp = sum( Transmit_power(:,1:t,:), 3);
figure()
plot( power_temp' );


figure()
power_temp = sum( Transmit_power(:,1:t,:), 3);
power_avg_all = bsxfun(@times, cumsum(power_temp, 2), 1./(1:t));
plot( mean(power_avg_all,1) )
plot( power_avg_all' );

%% Queues

queue_temp = Queues(:, 1:t) + Rate(:,1:t);
figure()
plot( queue_temp');
hold on
line([1 t], [threshold_Q threshold_Q]);
hold off

%% Queue samples
queue_tempx = maximum_queues(:, 1:t_sampling);
figure()
plot( queue_tempx');
hold on
line([1 t_sampling], [threshold_Q threshold_Q]);
line([1 t_sampling], [threshold_Q/2 threshold_Q/2]);
hold off


%% Power and queue
power_temp = sum( Transmit_power(:,1:t,:), 3);
queue_temp = Queues(:, 1:t) + Rate(:,1:t);
vq_reliable_temp = virtual_queue_reliability(:, 1:t);
vq_m1_temp = virtual_queue_moment_1(:, 1:t);
% vq_m2_temp = virtual_queue_moment_2(:,1:t);
figure()
for i=1:total_VUE_pairs
    subplot( 5, 1, 1);
    plot( power_temp(i,:) );
    title(['VUE pair = ' num2str(i) ....
        ' , with ' num2str(queue_exceed_counter(i)) ' extremes.']);
    ylabel('Power');
    subplot( 5, 1, 2);
    plot( queue_temp(i,:));
    ylabel('QSI');
    subplot( 5, 1, 3);
    plot( vq_reliable_temp(i,:));
    ylabel('VQ-reliable');
    subplot( 5, 1, 4);
    plot( vq_m1_temp(i,:));
    ylabel('VQ-M1');
    %     subplot( 5, 1, 5);
    %     plot( vq_m2_temp(i,:));
    %     ylabel('VQ-M2');
    pause(.5);
end

%% Power v Queue
power_temp = sum( Transmit_power(:,1:t,:), 3);
queue_temp = Queues(:, 1:t) + Rate(:,1:t);
figure()
plot( queue_temp(:), power_temp(:), '.');
xlabel('Queue');
ylabel('Power');

unreliability_vals = mean(Queues(:, 1:t) > threshold_Q, 2);
%%

stem(this_CINR_inv)
hold all
for i=1:length(water_level)
    plot(water_level(i)*ones(size(this_CINR_inv)));
end
plot(this_SI/tradeoff_lyapunov * ones(size(this_CINR_inv)), '-.')


%% FL
% valid_samples =[];
% for pair = 1:total_VUE_pairs
%     this_samplez = maximum_queues(pair, 1:t_sampling);
%     valid_samples = [valid_samples this_samplez( 1 + evtParam_global(3)*(this_samplez...
%         -evtParam_global(1))/evtParam_global(2) > 0 )];
%     plot(valid_samples, ones( size(valid_samples)), 'x');
% %     pause(1);
% end
valid_samples = sort(test_samples);
pd = fitdist( valid_samples, 'GeneralizedExtremeValue');
figure()
hold on
% observed distribution
[f_vals,x_vals]=ecdf(valid_samples);
myccdf = 1-f_vals;
plot(x_vals,myccdf, 'xb');
% matlab prediction
if (pd.k>=0)
    low_point_matlab = pd.mu - pd.sigma / pd.k;
    high_point_matlab = 10*floor(valid_samples(end)/10 + 1);
else
    low_point_matlab = 10*floor(valid_samples(1)/10 - 1);
    high_point_matlab = pd.mu - pd.sigma / pd.k;
end
xxx = linspace(low_point_matlab, high_point_matlab, 100);
yyy = 1 - cdf(pd, xxx);
plot(xxx,yyy,'-r');
% SGD approach
pdx = makedist('GeneralizedExtremeValue', 'mu', evtParam_global(1), ...
    'sigma', evtParam_global(2), 'k', evtParam_global(3) );
low_point_SGD = pdx.mu - pdx.sigma / pdx.k;
high_point_SGD = 10*floor(valid_samples(end)/10 + 1);
xx = linspace(low_point_SGD, high_point_SGD, 100);
yy = 1 - cdf(pdx, xxx);
plot(xx,yy,'-k');
hold off

%% Per user
for pair_val = 1:total_VUE_pairs
    these_samples = maximum_queues(pair_val, 1:t_sampling);
    pd = fitdist( these_samples', 'GeneralizedExtremeValue');
    disp(pd);
    pause(.1)
end

%% all users
these_samples = max( maximum_queues(:, 1:t_sampling), [], 1);
these_samples = sort(these_samples);
pd = fitdist( these_samples', 'GeneralizedExtremeValue');
disp(pd);
figure()
hold on
% observed distribution
[f_vals,x_vals]=ecdf(these_samples);
myccdf = 1-f_vals;
plot(x_vals,myccdf, 'xb');
% matlab prediction
if (pd.k>=0)
    low_point_matlab = pd.mu - pd.sigma / pd.k;
    high_point_matlab = 10*floor(these_samples(end)/10 + 1);
else
    low_point_matlab = 10*floor(these_samples(1)/10 - 1);
    high_point_matlab = pd.mu - pd.sigma / pd.k;
end
xxx = linspace(low_point_matlab, high_point_matlab, 100);
yyy = 1 - cdf(pd, xxx);
plot(xxx,yyy,'-r');
hold off


%% GPD

%% per user
for pair_val = 1:total_VUE_pairs
    these_samples_all = maximum_queues(pair_val, 1:t_sampling) - threshold_Q;
    these_samples = these_samples_all(these_samples_all>0);
    if length(these_samples)>2
        pd = fitdist( these_samples', 'generalizedpareto');
        disp(pair_val);
        disp(pd);
        pause(.1)
    end
end
%}
%% all
test_vals = {};%{'Method', 'Condition #1', 'Condition #1'};
these_samples_all = maximum_queues(:, 1:t_sampling) - threshold_Q;
these_samples = these_samples_all(these_samples_all>0);
these_samples = sort(these_samples);
if length(these_samples)>2
    pd = fitdist( these_samples, 'generalizedpareto');
    % disp('From Matlab:');
    % disp(pd);
    figure()
    subplot(4,2,[1 6]);
    hold on
    % observed distribution
    [f_vals,x_vals]=ecdf(these_samples);
    myccdf = 1-f_vals;
    plot(x_vals,myccdf, 'xb');
    legend_text = {'data points'};% by t = ' num2str(t)};
    % matlab prediction
    if (pd.k>=0)
        low_point_matlab = pd.theta - pd.sigma / pd.k;
        high_point_matlab = 10*floor(these_samples(end)/10 + 1);
    else
        low_point_matlab = 10*floor(these_samples(1)/10 - 1);
        high_point_matlab = pd.theta - pd.sigma / pd.k;
    end
    xxx = linspace(low_point_matlab, high_point_matlab, 100);
    yyy = 1 - cdf(pd, xxx);
    plot(xxx,yyy,'-r');
    legend_text{end+1} = ['Matlab: \sigma = ' num2str(pd.sigma) ', \xi = ' num2str(pd.k)];
    test_matlab = Test_MLE(pd.sigma,pd.k,these_samples);
    test_vals(end+1,:) = {'Matlab:', test_matlab(1), test_matlab(2)};
    % FL approach
    %{}
    pdx = makedist('generalizedpareto', 'theta', 0, ...
        'sigma', evtParam_global(1), 'k', evtParam_global(2) );
    % disp('From SGD:');
    % disp(pdx);
    if (pdx.k>=0)
        low_point_FL = pdx.theta - pdx.sigma / pdx.k;
        high_point_FL = 10*floor(these_samples(end)/10 + 1);
    else
        low_point_FL = 10*floor(these_samples(1)/10 - 1);
        high_point_FL = pdx.theta - pdx.sigma / pdx.k;
    end
    xx = linspace(low_point_FL, high_point_FL, 100);
    yy = 1 - cdf(pdx, xxx);
    plot(xx,yy,'-k');
    legend_text{end+1} = ['FL: \sigma = ' num2str(pdx.sigma) ', \xi = ' num2str(pdx.k)];
    test_FL = Test_MLE(pdx.sigma,pdx.k,these_samples);
    test_vals(end+1,:) = {'FL:', test_FL(1), test_FL(2)};
    %}
    %----------------------------
    %{
    gradient_global_org = [1 1000]; %dont make it zero
    evtParam_global_org = [50 0]; % [scale shape]
    local_sample_size = length(these_samples);
    gradient_SGD = zeros(size(gradient_global_org));
    count_VUE_samples = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % disp(pair);
    stepsize_SGD = stepsize_SVRG/local_sample_size;
    max_sample = max(these_samples);
    %             max_sample = max(local_sample_set);
    %             max_sample = local_sample_set;
    % Adjust parameters to fit into local data set
    % otherwise gradient becomes complex number
    evtParam_global_SGD = Projecting_Gradient_GPD( ...
        evtParam_global_org, max_sample );
    % localParam = evtParam_global;
    localParamSGD = evtParam_global_SGD - stepsize_SGD.* ...
        gradient_global_org;
    localParamSGD = Projecting_Gradient_GPD( ...
        localParamSGD, max_sample );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CHANGE
    randomized_samples = randperm(local_sample_size);
    % randomized_samples = 1:local_sample_size;
    for this_sample_index = randomized_samples
        this_sample = these_samples(this_sample_index);
        % disp(this_sample_index);
        % disp(this_sample*localParam(2) + localParam(1)>=0);
        % disp(this_sample*evtParam_global_this_set(2) + evtParam_global_this_set(1)>=0);
        gradient_now = ...
            Gradient_MaxLikelihood_GPD( localParamSGD, this_sample );
        gradient_now_global = ...
            Gradient_MaxLikelihood_GPD( evtParam_global_SGD, this_sample );
        localParam_unconstrained = localParamSGD - stepsize_SGD.*( ...
            gradient_now - gradient_now_global + gradient_global_org);
        % project to feasible set
        localParam_projected = Projecting_Gradient_GPD( ...
            localParam_unconstrained, max_sample );
        % If projection changes Param, gradient needs to be calculated.
        gradient_now = (localParamSGD - localParam_projected)./stepsize_SGD ...
            + gradient_now_global - gradient_global_org;
        localParamSGD = localParam_projected;
        %                         disp([this_sample_index this_sample localParam]);
        %                         disp(gradient_now);
        gradient_SGD = gradient_SGD + gradient_now;
        count_VUE_samples = count_VUE_samples + 1;
        %                         test_samples = [test_samples; this_sample];
    end %end for: randomize samples
    evtParam_local(pair,:) = localParamSGD;
    pdX = makedist('generalizedpareto', 'theta', 0, ...
        'sigma', localParamSGD(1), 'k', localParamSGD(2) );
    % disp('From SGD:');
    % disp(pdx);
    %{}
    if (pdX.k>=0)
        low_point_FL = pdX.theta - pdX.sigma / pdX.k;
        high_point_FL = 10*floor(these_samples(end)/10 + 1);
    else
        low_point_FL = 10*floor(these_samples(1)/10 - 1);
        high_point_FL = pdX.theta - pdX.sigma / pdX.k;
    end
    XX = linspace(low_point_FL, high_point_FL, 100);
    YY = 1 - cdf(pdX, xxx);
    plot(XX,YY,'-g');
    legend_text{end+1} = ['SGD-CEN: \sigma = ' num2str(pdX.sigma) ', \xi = ' num2str(pdX.k)];
    test_CEN = Test_MLE(pdX.sigma,pdX.k,these_samples);
    test_vals(end+1,:) = {'SGD-CEN:', test_CEN(1), test_CEN(2)};
    %}
    %----------------------------
    title( 'CCDF of excess values' )
    legend( legend_text );
    box on;
    hold off
    disp(test_vals)
    %{}
    Method = test_vals(:,1);
    Condition_1 = cell2mat(test_vals(:,2));
    Condition_2 = cell2mat(test_vals(:,3));
    T = table(Condition_1,Condition_2,'RowNames',Method);
    uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position', [0, 0, 1, .25]);
    %}
end