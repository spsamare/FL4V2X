% Network parameters

global total_VUE_pairs;
total_VUE_pairs = NetworkParams.TotalVUEPairs;%60; %multiplier of 4
global coherence_time bandwidth noise
coherence_time = 3 * 10^(-3); %%% 3 msec
bandwidth = 180 * 10^3;  %%% Bandwidth 180 kHz per RB
PSD = 10^(-174/10); %%% Power sepctral density -174 dBm/Hz
noise = PSD * bandwidth;


global total_resource_blocks max_transmit_power;
total_resource_blocks = NetworkParams.TotalResourceBlocks; % Number of RBs
max_transmit_power = (10^(10/10));%*total_resource_blocks;


global pathloss_coefficient pathloss_coefficient_NLOS pathloss_exponent
pathloss_exponent = 1.61;
pathloss_coefficient = 10^(-68.5/10);
pathloss_coefficient_NLOS = 10^(-54.5/10);
global distance_threshold_channel
distance_threshold_channel = 15; %meters, for WLOS



traffic_demand = NetworkParams.TrafficDemand*10^3; %%% 500 Kbps
average_packet_size = 500*8; %% bits per packet
average_packet_arrivals = traffic_demand/average_packet_size; %%packet/s