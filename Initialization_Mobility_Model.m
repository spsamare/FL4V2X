% Initialization_Mobility_Model
% Manhattan grid layout
% Change with caution.

%% non critical
offset_fixed = NetworkParams.DistanceFixed; %fixed TX-RX distance? if yes, use 1

global average_speed speed_fluctuation
average_speed = 60 * 1000/3600;                   %%% 60km/h
speed_fluctuation = 10 /100; %percentage of fluctuation over average speed

%% Critical: these values impact the mobility updates.
% %%% Coordinates of the x and y lanes are in 0, 120, and 240.
global offset_distance
offset_distance = NetworkParams.DistanceTxRx; %relies on the lane positions
