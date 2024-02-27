% Set_Parameters

% % % % % % % % % % % % % % % % % % %
% DelayModel:
%     0   Gamma
%     1   GammaTruncated
%     2   Uniform
% % % % % % % % % % % % % % % % % % %
NetworkParams = struct( ...
    'TotalVUEPairs', 80, ... % multiplier of 4
    'DistanceTxRx', 50, ... % in meters, on average 
    'DistanceFixed', 1, ... % gap between Tx and Rx is fixed?
    'TotalResourceBlocks', 60, ...
    'TrafficDemand', 500 ... %in kbps
    );

SimulationParams = struct( ...
    'SimulationTime', 30000,...
    'ClustersInclude', 1, ... %for cluster-based RB allocation
    'ClustersTotal', 10, ... 
    'ClustersTime', 100, ... % clusters update after this ms
    'IncludeReliability', 1, ... % 1 to yes 
    'InlcudeFederatedLearning', 1, ... % 1 to yes
    'QueueThresholdCoefficient', 1/8, ... % multiplied by 133(slot)*3ms*Arrival(bps)
    'QueueViolationProbability', 1e-3 ... % reliability condition
    );


ToolsParams = struct( ...
    'LyapunovTradeoff', 10^9, ... % highly depends on #VUEs
    'IterationsSVRGD', 50, ...% iteration within gradient decent
    'SamplingTime', 10, ...% time to sample max queue
    'FederatedLearningTime', 100 ... % FL takes after each of this time
    );
