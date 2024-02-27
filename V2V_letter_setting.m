%V2V_letter_setting
Avg_speed = 60 * 1000/3600;                   %%% 60km/h
acc_speed = 0;%10 /100; %percentage of fluctuation over average speed

PSD = 10^(-174/10);            %%% Power sepctral density -174 dBm/Hz
BW = 180 * 10^3;                %%% Bandwidth 180 kHz per RB
coherence_time = 3 * 10^(-3); %%% 3 msec
noise = PSD * BW;
%NumPair = 100;
NumGroup = 10;
%V = 0;
simutime = 1000000;
T0=100;
NumRB = 20;                        %%% Number of RBs
P_max = (10^(10/10))*NumRB;     %%% 10 dBm power budget
delta = 150; %%% Neighborhood region
zeta = 30;   %%% Similarity delay speed
triangle = 15;
Traffic_demand = 500*10^3; %%% 500 Kbps
avg_packet_size = 500*8;
BK_leng = BW * coherence_time;

A_thr_cen = 150*Traffic_demand/BW; %%% Mean network queue length bound,50(slot)*3ms*Arrival(bps)
B_thr_cen = 10^5; %%% Var network queue length bound
A_thr_loc = 150*Traffic_demand/BW; %%% Mean large local queue length bound
B_thr_loc = 10^5; %%% Var large local queue length bound,50(slot)*3ms*Arrival(bps)
%d = 100*Traffic_demand/BW; %%% for 2nd EVT theorem
viola_prob = 0.001; %%% Probabilistic constraint
Q_th =  100*Traffic_demand/BW; %%% Probabilistic constraint, 133(slot)*3ms*Arrival(bps)
POT_upper = 0.01; %%% Peak over threshold for shape parameter