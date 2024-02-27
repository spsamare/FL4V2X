function [...
    TX_location_next, ...
    RX_location_next, ...
    VUE_NS_direction_TX, ...
    VUE_EW_direction_TX, ...
    VUE_NS_direction_RX, ...
    VUE_EW_direction_RX ...
    ] = Update_VUE_Locations(...
    TX_location_current, ...
    RX_location_current, ...
    VUE_NS_direction_TX, ...
    VUE_EW_direction_TX, ...
    VUE_NS_direction_RX, ...
    VUE_EW_direction_RX, ...
    random_seed ...
    )

%% Globals
global average_speed speed_fluctuation coherence_time total_VUE_pairs


%% manupilate random state
current_random_state = rng;
rng(random_seed);

%% Updating location
instan_speed_TX = average_speed*...
    (1 + 2*(.5-rand(total_VUE_pairs,1))*speed_fluctuation);%%% meter/sec
instan_speed_RX = average_speed* ...
    (1 + 2*(.5-rand(total_VUE_pairs,1))*speed_fluctuation);%%% meter/sec
location_increase_TX = coherence_time * instan_speed_TX;
location_increase_RX = coherence_time * instan_speed_RX;


TX_location_next(:,1) =  TX_location_current(:,1) ...
    + VUE_EW_direction_TX.* location_increase_TX;
TX_location_next(:,2) =  TX_location_current(:,2) ....
    + VUE_NS_direction_TX.* location_increase_TX;
RX_location_next(:,1) =  RX_location_current(:,1) ...
    + VUE_EW_direction_RX.* location_increase_RX;
RX_location_next(:,2) =  RX_location_current(:,2) ...
    + VUE_NS_direction_RX.* location_increase_RX;

%% Updating while near intersection for TX
dummy_ind_1 = (find( ( (TX_location_next(:,1)    ) .* ...
    (TX_location_current(:,1)    ) ) < 0 )).';
dummy_ind_2 = (find( ( (TX_location_next(:,1)-120) .* ....
    (TX_location_current(:,1)-120) ) < 0 )).';
dummy_ind_3 = (find( ( (TX_location_next(:,1)-240) .* ...
    (TX_location_current(:,1)-240) ) < 0 )).';
dummy_ind_4 = (find( ( (TX_location_next(:,2)    ) .* ....
    (TX_location_current(:,2)    ) ) < 0 )).';
dummy_ind_5 = (find( ( (TX_location_next(:,2)-120) .* ....
    (TX_location_current(:,2)-120) ) < 0 )).';
dummy_ind_6 = (find( ( (TX_location_next(:,2)-240) .* ...
    (TX_location_current(:,2)-240) ) < 0 )).';

for k = dummy_ind_1
    if TX_location_current(k,2)==0
        TX_location_next(k,1) = 0;
        TX_location_next(k,2) = location_increase_TX(k) ...
            - TX_location_current(k,1);
        VUE_EW_direction_TX(k,1) = 0;
        VUE_NS_direction_TX(k,1) = 1;
    elseif TX_location_current(k,2)== 120
        TX_location_next(k,1) = 0;
        dummy_direction = 2*round(rand(1))-1;
        TX_location_next(k,2) = 120 + dummy_direction*(location_increase_TX(k)...
            - TX_location_current(k,1));
        VUE_EW_direction_TX(k,1) = 0;
        VUE_NS_direction_TX(k,1) = dummy_direction;
    elseif TX_location_current(k,2)== 240
        TX_location_next(k,1) = 0;
        TX_location_next(k,2) = 240 - (location_increase_TX(k)...
            - TX_location_current(k,1));
        VUE_EW_direction_TX(k,1) = 0;
        VUE_NS_direction_TX(k,1) = -1;
    end
end


for k = dummy_ind_2
    if (TX_location_current(k,1)< 120)
        if TX_location_current(k,2)==0
            dummy_direction = round(rand(1));
            VUE_EW_direction_TX(k,1) = mod(VUE_EW_direction_TX(k,1)...
                +dummy_direction,2);
            VUE_NS_direction_TX(k,1) = dummy_direction;
            reman = location_increase_TX(k) + TX_location_current(k,1) - 120;
            TX_location_next(k,1) = 120 + VUE_EW_direction_TX(k,1)*reman;
            TX_location_next(k,2) =  VUE_NS_direction_TX(k,1)*reman;
            
        elseif TX_location_current(k,2)== 120
            dummy_direction = ceil(3*rand(1))-2;
            if (dummy_direction==1) || (dummy_direction == -1)
                VUE_EW_direction_TX(k,1) = 0;
                VUE_NS_direction_TX(k,1) = dummy_direction;
            end
            reman = location_increase_TX(k) + TX_location_current(k,1)-120;
            TX_location_next(k,1) = 120 + VUE_EW_direction_TX(k,1)*reman;
            TX_location_next(k,2) = 120 + VUE_NS_direction_TX(k,1)*reman;
            
        elseif TX_location_current(k,2) == 240
            dummy_direction = round(rand(1))-1;
            VUE_EW_direction_TX(k,1) = mod(VUE_EW_direction_TX(k,1)...
                +dummy_direction,2);
            VUE_NS_direction_TX(k,1) = dummy_direction;
            reman = location_increase_TX(k) + TX_location_current(k,1)-120;
            TX_location_next(k,1) = 120 + VUE_EW_direction_TX(k,1)*reman;
            TX_location_next(k,2) =  240 + VUE_NS_direction_TX(k,1)*reman;
        end
    elseif (TX_location_current(k,1)> 120)
        if TX_location_current(k,2)==0
            dummy_direction = round(rand(1));
            VUE_EW_direction_TX(k,1) = VUE_EW_direction_TX(k,1)...
                +dummy_direction;
            VUE_NS_direction_TX(k,1) = dummy_direction;
            reman = location_increase_TX(k) - (TX_location_current(k,1)-120);
            TX_location_next(k,1) = 120 + VUE_EW_direction_TX(k,1)*reman;
            TX_location_next(k,2) =  VUE_NS_direction_TX(k,1)*reman;
            
        elseif TX_location_current(k,2)== 120
            dummy_direction = ceil(3*rand(1))-2;
            if (dummy_direction==1) || (dummy_direction == -1)
                VUE_EW_direction_TX(k,1) = 0;
                VUE_NS_direction_TX(k,1) = dummy_direction;
            end
            reman = location_increase_TX(k) - TX_location_current(k,1)+120;
            TX_location_next(k,1) = 120 + VUE_EW_direction_TX(k,1)*reman;
            TX_location_next(k,2) = 120 + VUE_NS_direction_TX(k,1)*reman;
            
        elseif TX_location_current(k,2)== 240
            dummy_direction = round(rand(1))-1;
            VUE_EW_direction_TX(k,1) = mod(VUE_EW_direction_TX(k,1)...
                +dummy_direction,-2);
            VUE_NS_direction_TX(k,1) = dummy_direction;
            reman = location_increase_TX(k) - TX_location_current(k,1)+120;
            TX_location_next(k,1) = 120 + VUE_EW_direction_TX(k,1)*reman;
            TX_location_next(k,2) = 240 + VUE_NS_direction_TX(k,1)*reman;
        end
    end
end


for k = dummy_ind_3
    if TX_location_current(k,2)==0
        TX_location_next(k,1) = 240;
        TX_location_next(k,2) = location_increase_TX(k) ...
            - (240-TX_location_current(k,1));
        VUE_EW_direction_TX(k,1) = 0;
        VUE_NS_direction_TX(k,1) = 1;
    elseif TX_location_current(k,2)== 120
        TX_location_next(k,1) = 240;
        dummy_direction = 2*round(rand(1))-1;
        TX_location_next(k,2) = 120 + dummy_direction*(location_increase_TX(k)...
            - (240-TX_location_current(k,1)) );
        VUE_EW_direction_TX(k,1) = 0;
        VUE_NS_direction_TX(k,1) = dummy_direction;
    elseif TX_location_current(k,2)== 240
        TX_location_next(k,1) = 240;
        TX_location_next(k,2) = 240 - (location_increase_TX(k)...
            - (240-TX_location_current(k,1)));
        VUE_EW_direction_TX(k,1) = 0;
        VUE_NS_direction_TX(k,1) = -1;
    end
end





for k = dummy_ind_4
    if TX_location_current(k,1)==0
        TX_location_next(k,2) = 0;
        TX_location_next(k,1) = location_increase_TX(k)...
            - TX_location_current(k,2);
        VUE_EW_direction_TX(k,1) = 1;
        VUE_NS_direction_TX(k,1) = 0;
    elseif TX_location_current(k,1)== 120
        TX_location_next(k,2) = 0;
        dummy_direction = 2*round(rand(1))-1;
        TX_location_next(k,1) = 120 + dummy_direction*(location_increase_TX(k)...
            - TX_location_current(k,2));
        VUE_EW_direction_TX(k,1) = dummy_direction;
        VUE_NS_direction_TX(k,1) = 0;
    elseif TX_location_current(k,1)== 240
        TX_location_next(k,2) = 0;
        TX_location_next(k,1) = 240 - (location_increase_TX(k)...
            - TX_location_current(k,2));
        VUE_EW_direction_TX(k,1) = -1;
        VUE_NS_direction_TX(k,1) = 0;
    end
end


for k = dummy_ind_5
    if (TX_location_current(k,2)< 120)
        if TX_location_current(k,1)==0
            dummy_direction = round(rand(1));
            VUE_NS_direction_TX(k,1) = mod(VUE_NS_direction_TX(k,1)...
                +dummy_direction,2);
            VUE_EW_direction_TX(k,1) = dummy_direction;
            reman = location_increase_TX(k) + TX_location_current(k,2) - 120;
            TX_location_next(k,1) =  VUE_EW_direction_TX(k,1)*reman;
            TX_location_next(k,2) = 120 + VUE_NS_direction_TX(k,1)*reman;
            
        elseif TX_location_current(k,1)== 120
            dummy_direction = ceil(3*rand(1))-2;
            if (dummy_direction==1) || (dummy_direction == -1)
                VUE_NS_direction_TX(k,1) = 0;
                VUE_EW_direction_TX(k,1) = dummy_direction;
            end
            reman = location_increase_TX(k) + TX_location_current(k,2)-120;
            TX_location_next(k,1) = 120 + VUE_EW_direction_TX(k,1)*reman;
            TX_location_next(k,2) = 120 + VUE_NS_direction_TX(k,1)*reman;
            
        elseif TX_location_current(k,1)== 240
            dummy_direction = round(rand(1))-1;
            VUE_NS_direction_TX(k,1) = mod(VUE_NS_direction_TX(k,1)...
                +dummy_direction,2);
            VUE_EW_direction_TX(k,1) = dummy_direction;
            reman = location_increase_TX(k) + TX_location_current(k,2)-120;
            TX_location_next(k,1) = 240 + VUE_EW_direction_TX(k,1)*reman;
            TX_location_next(k,2) =  120 + VUE_NS_direction_TX(k,1)*reman;
        end
        
    elseif (TX_location_current(k,2)> 120)
        
        if TX_location_current(k,1)==0
            dummy_direction = round(rand(1));
            VUE_NS_direction_TX(k,1) = VUE_NS_direction_TX(k,1)...
                +dummy_direction;
            VUE_EW_direction_TX(k,1) = dummy_direction;
            reman = location_increase_TX(k) - (TX_location_current(k,2)-120);
            TX_location_next(k,1) =  VUE_EW_direction_TX(k,1)*reman;
            TX_location_next(k,2) = 120 + VUE_NS_direction_TX(k,1)*reman;
            
        elseif TX_location_current(k,1)== 120
            dummy_direction = ceil(3*rand(1))-2;
            if (dummy_direction==1) || (dummy_direction == -1)
                VUE_NS_direction_TX(k,1) = 0;
                VUE_EW_direction_TX(k,1) = dummy_direction;
            end
            reman = location_increase_TX(k) - TX_location_current(k,2)+120;
            TX_location_next(k,1) = 120 + VUE_EW_direction_TX(k,1)*reman;
            TX_location_next(k,2) = 120 + VUE_NS_direction_TX(k,1)*reman;
            
        elseif TX_location_current(k,1)== 240
            dummy_direction = round(rand(1))-1;
            VUE_NS_direction_TX(k,1) = mod(VUE_NS_direction_TX(k,1)...
                +dummy_direction,-2);
            VUE_EW_direction_TX(k,1) = dummy_direction;
            reman = location_increase_TX(k) - TX_location_current(k,2)+120;
            TX_location_next(k,1) = 240 + VUE_EW_direction_TX(k,1)*reman;
            TX_location_next(k,2) =  120 + VUE_NS_direction_TX(k,1)*reman;
        end
    end
end



for k = dummy_ind_6
    if TX_location_current(k,1)==0
        TX_location_next(k,2) = 240;
        TX_location_next(k,1) = location_increase_TX(k)...
            - (240-TX_location_current(k,2));
        VUE_EW_direction_TX(k,1) = 1;
        VUE_NS_direction_TX(k,1) = 0;
    elseif TX_location_current(k,1)== 120
        TX_location_next(k,2) = 240;
        dummy_direction = 2*round(rand(1))-1;
        TX_location_next(k,1) = 120 + dummy_direction*(location_increase_TX(k)...
            - (240-TX_location_current(k,2)) );
        VUE_NS_direction_TX(k,1) = 0;
        VUE_EW_direction_TX(k,1) = dummy_direction;
    elseif TX_location_current(k,1)== 240
        TX_location_next(k,2) = 240;
        TX_location_next(k,1) = 240 - (location_increase_TX(k)...
            - (240-TX_location_current(k,2)) );
        VUE_EW_direction_TX(k,1) = -1;
        VUE_NS_direction_TX(k,1) = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Updating while near intersection for RX
dummy_ind_1 = (find( ( (RX_location_next(:,1)    ) .*...
    (RX_location_current(:,1)    ) ) < 0 )).';
dummy_ind_2 = (find( ( (RX_location_next(:,1)-120) .*...
    (RX_location_current(:,1)-120) ) < 0 )).';
dummy_ind_3 = (find( ( (RX_location_next(:,1)-240) .*...
    (RX_location_current(:,1)-240) ) < 0 )).';
dummy_ind_4 = (find( ( (RX_location_next(:,2)    ) .*...
    (RX_location_current(:,2)    ) ) < 0 )).';
dummy_ind_5 = (find( ( (RX_location_next(:,2)-120) .*...
    (RX_location_current(:,2)-120) ) < 0 )).';
dummy_ind_6 = (find( ( (RX_location_next(:,2)-240) .*...
    (RX_location_current(:,2)-240) ) < 0 )).';



for k = dummy_ind_1
    if RX_location_current(k,2)==0
        RX_location_next(k,1) = 0;
        RX_location_next(k,2) = location_increase_RX(k) - RX_location_current(k,1);
        VUE_EW_direction_RX(k,1) = 0;
        VUE_NS_direction_RX(k,1) = 1;
    elseif RX_location_current(k,2)== 120
        RX_location_next(k,1) = 0;
        RX_location_next(k,2) = 120 + VUE_NS_direction_TX(k,1)*(location_increase_RX(k)...
            - RX_location_current(k,1));
        VUE_EW_direction_RX(k,1) = 0;
        VUE_NS_direction_RX(k,1) = VUE_NS_direction_TX(k,1);
        
    elseif RX_location_current(k,2)== 240
        RX_location_next(k,1) = 0;
        RX_location_next(k,2) = 240 - (location_increase_RX(k)...
            - RX_location_current(k,1));
        VUE_EW_direction_RX(k,1) = 0;
        VUE_NS_direction_RX(k,1) = -1;
    end
end



for k = dummy_ind_2
    if (RX_location_current(k,1)< 120)
        if RX_location_current(k,2)==0
            VUE_EW_direction_RX(k,1) = VUE_EW_direction_TX(k,1);
            VUE_NS_direction_RX(k,1) = VUE_NS_direction_TX(k,1);
            reman = location_increase_RX(k) + RX_location_current(k,1) - 120;
            RX_location_next(k,1) = 120 + VUE_EW_direction_RX(k,1)*reman;
            RX_location_next(k,2) =  VUE_NS_direction_RX(k,1)*reman;
            
        elseif RX_location_current(k,2)== 120
            
            VUE_EW_direction_RX(k,1) = VUE_EW_direction_TX(k,1);
            VUE_NS_direction_RX(k,1) = VUE_NS_direction_TX(k,1);
            reman = location_increase_RX(k) + RX_location_current(k,1)-120;
            RX_location_next(k,1) = 120 + VUE_EW_direction_RX(k,1)*reman;
            RX_location_next(k,2) = 120 + VUE_NS_direction_RX(k,1)*reman;
            
        elseif RX_location_current(k,2)== 240
            
            VUE_EW_direction_RX(k,1) = VUE_EW_direction_TX(k,1);
            VUE_NS_direction_RX(k,1) = VUE_NS_direction_TX(k,1);
            reman = location_increase_RX(k) + RX_location_current(k,1)-120;
            RX_location_next(k,1) = 120 + VUE_EW_direction_RX(k,1)*reman;
            RX_location_next(k,2) =  240 + VUE_NS_direction_RX(k,1)*reman;
        end
    elseif (RX_location_current(k,1)> 120)
        if RX_location_current(k,2)==0
            VUE_EW_direction_RX(k,1) = VUE_EW_direction_TX(k,1);
            VUE_NS_direction_RX(k,1) = VUE_NS_direction_TX(k,1);
            reman = location_increase_RX(k) - (RX_location_current(k,1)-120);
            RX_location_next(k,1) = 120 + VUE_EW_direction_RX(k,1)*reman;
            RX_location_next(k,2) =  VUE_NS_direction_RX(k,1)*reman;
            
        elseif RX_location_current(k,2)== 120
            VUE_EW_direction_RX(k,1) = VUE_EW_direction_TX(k,1);
            VUE_NS_direction_RX(k,1) = VUE_NS_direction_TX(k,1);
            reman = location_increase_RX(k) - RX_location_current(k,1)+120;
            RX_location_next(k,1) = 120 + VUE_EW_direction_RX(k,1)*reman;
            RX_location_next(k,2) = 120 + VUE_NS_direction_RX(k,1)*reman;
            
        elseif RX_location_current(k,2)== 240
            VUE_EW_direction_RX(k,1) = VUE_EW_direction_TX(k,1);
            VUE_NS_direction_RX(k,1) = VUE_NS_direction_TX(k,1);
            reman = location_increase_RX(k) - RX_location_current(k,1)+120;
            RX_location_next(k,1) = 120 + VUE_EW_direction_RX(k,1)*reman;
            RX_location_next(k,2) =  240 + VUE_NS_direction_RX(k,1)*reman;
        end
    end
end




for k = dummy_ind_3
    if RX_location_current(k,2)==0
        RX_location_next(k,1) = 240;
        RX_location_next(k,2) = location_increase_RX(k)...
            - (240-RX_location_current(k,1));
        VUE_EW_direction_RX(k,1) = 0;
        VUE_NS_direction_RX(k,1) = 1;
    elseif RX_location_current(k,2)== 120
        RX_location_next(k,1) = 240;
        VUE_EW_direction_RX(k,1) = 0;
        VUE_NS_direction_RX(k,1) = VUE_NS_direction_TX(k,1);
        RX_location_next(k,2) = 120 + VUE_NS_direction_RX(k,1)*(location_increase_RX(k)...
            - (240-RX_location_current(k,1)) );
        
    elseif RX_location_current(k,2)== 240
        RX_location_next(k,1) = 240;
        RX_location_next(k,2) = 240 - (location_increase_RX(k)...
            - (240-RX_location_current(k,1)));
        VUE_EW_direction_RX(k,1) = 0;
        VUE_NS_direction_RX(k,1) = -1;
    end
end





for k = dummy_ind_4
    if RX_location_current(k,1)==0
        RX_location_next(k,2) = 0;
        RX_location_next(k,1) = location_increase_RX(k)...
            - RX_location_current(k,2);
        VUE_EW_direction_RX(k,1) = 1;
        VUE_NS_direction_RX(k,1) = 0;
    elseif RX_location_current(k,1)== 120
        VUE_EW_direction_RX(k,1) = VUE_EW_direction_TX(k,1);
        VUE_NS_direction_RX(k,1) = 0;
        RX_location_next(k,2) = 0;
        RX_location_next(k,1) = 120 + VUE_EW_direction_RX(k,1) *(location_increase_RX(k)...
            - RX_location_current(k,2));
        
    elseif RX_location_current(k,1)== 240
        RX_location_next(k,2) = 0;
        RX_location_next(k,1) = 240 - (location_increase_RX(k)...
            - RX_location_current(k,2));
        VUE_EW_direction_RX(k,1) = -1;
        VUE_NS_direction_RX(k,1) = 0;
    end
end


for k = dummy_ind_5
    if RX_location_current(k,2)< 120
        
        if RX_location_current(k,1)==0
            VUE_NS_direction_RX(k,1) = VUE_NS_direction_TX(k,1);
            VUE_EW_direction_RX(k,1) = VUE_EW_direction_TX(k,1);
            reman = location_increase_RX(k) + RX_location_current(k,2) - 120;
            RX_location_next(k,1) =  VUE_EW_direction_RX(k,1)*reman;
            RX_location_next(k,2) = 120 + VUE_NS_direction_RX(k,1)*reman;
        elseif RX_location_current(k,1)== 120
            VUE_NS_direction_RX(k,1) = VUE_NS_direction_TX(k,1);
            VUE_EW_direction_RX(k,1) = VUE_EW_direction_TX(k,1);
            reman = location_increase_RX(k) + RX_location_current(k,2)-120;
            RX_location_next(k,1) = 120 + VUE_EW_direction_RX(k,1)*reman;
            RX_location_next(k,2) = 120 + VUE_NS_direction_RX(k,1)*reman;
            
        elseif RX_location_current(k,1)== 240
            VUE_NS_direction_RX(k,1) = VUE_NS_direction_TX(k,1);
            VUE_EW_direction_RX(k,1) = VUE_EW_direction_TX(k,1);
            reman = location_increase_RX(k) + RX_location_current(k,2)-120;
            RX_location_next(k,1) = 240 + VUE_EW_direction_RX(k,1)*reman;
            RX_location_next(k,2) =  120 + VUE_NS_direction_RX(k,1)*reman;
        end
        
    elseif RX_location_current(k,2)> 120
        if RX_location_current(k,1)==0
            VUE_NS_direction_RX(k,1) = VUE_NS_direction_TX(k,1);
            VUE_EW_direction_RX(k,1) = VUE_EW_direction_TX(k,1);
            reman = location_increase_RX(k) - (RX_location_current(k,2)-120);
            RX_location_next(k,1) =  VUE_EW_direction_RX(k,1)*reman;
            RX_location_next(k,2) = 120 + VUE_NS_direction_RX(k,1)*reman;
            
        elseif RX_location_current(k,1)== 120
            VUE_NS_direction_RX(k,1) = VUE_NS_direction_TX(k,1);
            VUE_EW_direction_RX(k,1) = VUE_EW_direction_TX(k,1);
            reman = location_increase_RX(k) - RX_location_current(k,2)+120;
            RX_location_next(k,1) = 120 + VUE_EW_direction_RX(k,1)*reman;
            RX_location_next(k,2) = 120 + VUE_NS_direction_RX(k,1)*reman;
            
        elseif RX_location_current(k,1)== 240
            VUE_NS_direction_RX(k,1) = VUE_NS_direction_TX(k,1) ;
            VUE_EW_direction_RX(k,1) =  VUE_EW_direction_TX(k,1);
            reman = location_increase_RX(k) - RX_location_current(k,2)+120;
            RX_location_next(k,1) = 240 + VUE_EW_direction_RX(k,1)*reman;
            RX_location_next(k,2) =  120 + VUE_NS_direction_RX(k,1)*reman;
        end
    end
end



for k = dummy_ind_6
    if RX_location_current(k,1)==0
        RX_location_next(k,2) = 240;
        RX_location_next(k,1) = location_increase_RX(k)...
            - (240-RX_location_current(k,2));
        VUE_EW_direction_RX(k,1) = 1;
        VUE_NS_direction_RX(k,1) = 0;
    elseif RX_location_current(k,1)== 120
        RX_location_next(k,2) = 240;
        VUE_NS_direction_RX(k,1) = VUE_NS_direction_TX(k,1);
        VUE_EW_direction_RX(k,1) = VUE_EW_direction_TX(k,1);
        RX_location_next(k,1) = 120 + VUE_EW_direction_RX(k,1)*(location_increase_RX(k)...
            - (240-RX_location_current(k,2)) );
        
    elseif RX_location_current(k,1)== 240
        RX_location_next(k,2) = 240;
        RX_location_next(k,1) = 240 - (location_increase_RX(k) -....
            (240-RX_location_current(k,2)));
        VUE_EW_direction_RX(k,1) = -1;
        VUE_NS_direction_RX(k,1) = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% restore original random state
rng(current_random_state);

end