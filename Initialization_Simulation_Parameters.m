rng(RANDOM_SEED);
simulation_time = SimulationParams.SimulationTime;
random_seeds_location = randperm(simulation_time);
random_seeds_channel = randperm(simulation_time);

show_figures = 0; %yes = 1, no = 0
show_directions = 0; %yes = 1, no = 0


% switching between baseline models and proposed model
include_reliability_constraint = SimulationParams.IncludeReliability; %yes = 1, no = 0
include_federated_learning = SimulationParams.InlcudeFederatedLearning; %yes = 1, no = 0




% store arrivals for each VUE pair over time
Arrivals = zeros( total_VUE_pairs, simulation_time);
% Transmit_power = zeros( total_VUE_pairs, simulation_time, total_resource_blocks); % issue of saving large matrix
Transmit_power = zeros( total_VUE_pairs, simulation_time); 
Capacity = zeros( total_VUE_pairs, simulation_time, total_resource_blocks);
Rate = zeros( total_VUE_pairs, simulation_time); % depends on the queue as well
Queues = zeros( total_VUE_pairs, simulation_time);
% not from scratch
% Queues(:,1) = exprnd( 50*traffic_demand/bandwidth, [total_VUE_pairs, 1]);

% -------------------------------------
% Clustering parameter
perform_cluster_based_RBs = SimulationParams.ClustersInclude;
perform_random_RBs = 1 * (perform_cluster_based_RBs==0);

resource_allocation_table = ones( total_VUE_pairs, total_resource_blocks);
t_clustering = SimulationParams.ClustersTime;
global total_clusters neighborhood similarity_std
total_clusters = SimulationParams.ClustersTotal;
neighborhood = 150; %%% Neighborhood region [m]
similarity_std = 30/sqrt(2);   %%% standard deviation of Similarity

% Lyapunov parameters
virtual_queue_reliability = zeros( total_VUE_pairs, simulation_time);
virtual_queue_moment_1 = zeros( total_VUE_pairs, simulation_time);
%------removed for POT
% virtual_queue_moment_2 = zeros( total_VUE_pairs, simulation_time);

violation_probability = SimulationParams.QueueViolationProbability; %%% Probabilistic constraint
threshold_Q =  SimulationParams.QueueThresholdCoefficient *...
    100*traffic_demand/bandwidth; %%% Probabilistic constraint, 133(slot)*3ms*Arrival(bps)
global tradeoff_lyapunov
tradeoff_lyapunov = ToolsParams.LyapunovTradeoff; %%%test 10^7

% -------------------------------------
% EVT and FL parameters
counter_federated = 0;

stepsize_SVRG = [50 10]; % step size of gradient decent
iterations_SVGR = ToolsParams.IterationsSVRGD;
max_sample_global = 1;

window_sampling = ToolsParams.SamplingTime; % pick a smaple within a window
window_federated = ToolsParams.FederatedLearningTime; % share models after this period (in terms of samples)

t_sampling = 0;
% t_federated = 0; %not used

maximum_queues = zeros( total_VUE_pairs, ceil(simulation_time/window_sampling));
%global models
% gradient_global = ones(1, 3); %dont make it zero
% evtParam_global = [250 2 1]; % [location scale shape]
% modified for GPD
gradient_global = [1 1000]; %dont make it zero
evtParam_global = [50 0]; % [scale shape]
samples_global = 0;
%local models
gradient_local = repmat( gradient_global, [total_VUE_pairs, 1]); %
evtParam_local = repmat( evtParam_global, [total_VUE_pairs, 1]); % [location scale shape]
samples_local = zeros(total_VUE_pairs, 1);
% store during local updates only
% active_learning_indicator = ones(total_VUE_pairs, 1);
% gradient_VUE_dummy = zeros(size(gradient_local));
% count_samples = zeros(total_VUE_pairs, 1);

threshold_maxQ_moment_1 = threshold_Q + ...
    evtParam_global(1) / (1 - evtParam_global(2)); %%% Peak over threshold
%---------changed to POT-----------
% threshold_maxQ_moment_1 = 150*traffic_demand/bandwidth; %%% Mean large local queue length bound
% threshold_maxQ_moment_2 = 10^5; %%% Var large local queue length bound,50(slot)*3ms*Arrival(bps)

queue_exceed_counter = zeros( total_VUE_pairs, 1);


% interference estimation
interferenceNnoice_estimation = noise * ...
    ones( total_VUE_pairs, total_resource_blocks);
weight_old_est = .5;