close all

% observed distribution
%{
[f_vals,x_vals]=ecdf(these_samples);
myccdf = 1-f_vals;
plot(x_vals,myccdf, 'xb');
title( ['Matlab: \sigma = ' num2str(pd.sigma) ', \xi = ' num2str(pd.k)] );
legend_text = {'data points'};
legend( legend_text );
box on;
%}
count = 1;

temp = []; tempGrad = []; 

%----------------------------
gradient_global_org = [1 1000]; %dont make it zero
evtParam_global_org = [50 0]; % [scale shape]
local_sample_size = length(these_samples);
gradient_SGD = zeros(size(gradient_global_org));
count_VUE_samples = 0;
tempParam = evtParam_global_org;
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

myCon = 1;
gradient_global_old = gradient_global_org;
% gradient_SGD = [0 0];

while( count<500 )
% while( myCon )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CHANGE
    gradient_SGD = [0 0];
    randomized_samples = randperm(local_sample_size);
    % randomized_samples = 1:local_sample_size;
    for this_sample_index = randomized_samples
        this_sample = these_samples(this_sample_index);
        gradient_now = ...
            Gradient_MaxLikelihood_GPD( localParamSGD, this_sample );
        gradient_now_global = ...
            Gradient_MaxLikelihood_GPD( evtParam_global_SGD, this_sample );
        localParam_unconstrained = localParamSGD - stepsize_SGD.*( ...
            gradient_now - gradient_now_global + gradient_global_old);
        % project to feasible set
        localParam_projected = Projecting_Gradient_GPD( ...
            localParam_unconstrained, max_sample );
        % If projection changes Param, gradient needs to be calculated.
%         gradient_now = (localParamSGD - localParam_projected)./stepsize_SGD ...
%             + gradient_now_global - gradient_global_old;
        localParamSGD = localParam_projected;
        gradient_SGD = gradient_SGD + gradient_now;
% % % % % % % % %         gradient_SGD = gradient_SGD + gradient_now_global;
        count_VUE_samples = count_VUE_samples + 1;
    end %end for: randomize samples
%     evtParam_local(pair,:) = localParamSGD;
    gradient_global_old = gradient_SGD/local_sample_size;
%     gradient_global_old = gradient_SGD/count_VUE_samples;
    evtParam_global_SGD = localParamSGD;
    
    pdX = makedist('generalizedpareto', 'theta', 0, ...
        'sigma', localParamSGD(1), 'k', localParamSGD(2) );
    %{
    if (pdX.k>=0)
        low_point_SGD = pdX.theta - pdX.sigma / pdX.k;
        high_point_SGD = 10*floor(these_samples(end)/10 + 1);
    else
        low_point_SGD = 10*floor(these_samples(1)/10 - 1);
        high_point_SGD = pdX.theta - pdX.sigma / pdX.k;
    end
    XX = linspace(low_point_SGD, high_point_SGD, 100);
    YY = 1 - cdf(pdX, xxx);
    plot(XX,YY,'-g');
    %----------------------------
    
    legend_text{end+1} = ['SGD-CEN: \sigma = ' num2str(pdX.sigma) ', \xi = ' num2str(pdX.k)];
    legend(legend_text);
    
    disp('Initial gradient: ');
    disp(gradient_global_org);
    disp('Final gradient: ');
    disp(gradient_SGD.*stepsize_SGD);
    %}
    temp = [temp; gradient_SGD.*stepsize_SGD];
    tempGrad = [tempGrad; gradient_SGD];
    tempParam = [tempParam; localParamSGD];
    
    
%     plot(count, pd.sigma, 'ro');
%     plot(count, -10*pd.k, 'bo');
%     plot(count, pdX.sigma, 'k.');
%     plot(count, -10*pdX.k, 'g.');
% plot(count, pd.sigma - pdX.sigma, 'rx');

    b = pdX.sigma; a = pdX.k;
    val11 = sum( ( these_samples / b )./( 1 + a * these_samples / b ) ) - length(these_samples)/ (1 + a);
    val22 = sum( log( 1 + a * these_samples / b ) ) - length(these_samples)*a;
    
    figure(1)
    hold on
    plot( count, val11, 'rx');
    
    figure(2)
    hold on
    plot( count, val22, 'bx');
    
    
    count = count + 1;
    pause(0.01);
    
%     x_inp = input('For another round, press 1 and enter: ');
%     if(x_inp~=1)
%         myCon = 0;
%     end
    
end
% ---------------
hold off







%%
ARRIVALS = Arrivals(:,1:t);
save( 'testMod0.mat', ...
    'ARRIVALS', ...
    'RX_location_center', ...
    'TX_location_center', ...
    'VUE_EW_direction_RX', ...
    'VUE_EW_direction_TX', ...
    'VUE_NS_direction_RX', ...
    'VUE_NS_direction_TX' ...
    );



