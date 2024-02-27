% Mod_x
% x = 0 : no reliability, no EVT
% x = 1 : no EVT
% x = 2 : proposed

myModel = include_reliability_constraint + include_federated_learning;

filename = [...
    'results_MOD_' num2str(myModel) '_' ...
    'VUE_' num2str(total_VUE_pairs, '%04d') '_'...
    'Top_' num2str(RANDOM_SEED, '%03d') ...
    '.mat'...
    ];

save( filename, ...
    't', ...
    'total_VUE_pairs', ...
    'offset_distance', ...
    'traffic_demand', ...
    'threshold_Q', ...
    'include_reliability_constraint', ...
    'include_federated_learning', ...
    'violation_probability', ...
    'window_sampling', ...
    'window_federated', ...
    'tradeoff_lyapunov', ...
    'Queues', ...
    'Rate', ...
    'Transmit_power', ...
    'maximum_queues', ...
    'evtParam_global', ...
    'counter_federated' ...
    );
disp(['Data saving completed at time = ' num2str(t)]); %%