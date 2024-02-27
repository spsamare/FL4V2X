% function V2V_EVT(NumPair,V) %140 pairs at max
NumPair = 16;
V = 0;


V2V_letter_setting;
% simutime=10000;
%%%% Coordinates of the x and y lanes are in 0, 120, and 240.
Mean_Pow = zeros(NumPair,3);

% Queue_cen = zeros(NumPair,1);
% NetQeue_A = 0;
% NetQeue_B = 0;
% Mean_NetQeue_A = zeros(1,simutime+1);
% Mean_NetQeue_B = zeros(1,simutime+1);
% Rate_cen = zeros(NumPair,simutime);
% Actual_rate_cen = zeros(NumPair,simutime);
    
Queue_loc = zeros(NumPair,1);
LocQeue_A = zeros(NumPair,1);
LocQeue_B = zeros(NumPair,1);
Mean_LocQeue_A = zeros(NumPair,simutime+1);
Mean_LocQeue_B = zeros(NumPair,simutime+1);
mu = zeros(NumPair,1);
sigma = zeros(NumPair,1);
xi = zeros(NumPair,1);
Rate_loc = zeros(NumPair,simutime);
Actual_rate_loc = zeros(NumPair,simutime);

Queue_base = zeros(NumPair,1);
BaseQeue_A = zeros(NumPair,1);
Mean_BaseQeue_A = zeros(NumPair,simutime+1);
Rate_base = zeros(NumPair,simutime);
Actual_rate_base = zeros(NumPair,simutime);

Arrival = zeros(NumPair,simutime);
%
%
%%%%%%%%%%%%%%%%%%%%%%% Initializing the vehicles' positions

TX_location_curr(1:NumPair/4,1) = 120*ceil(3*rand(NumPair/4,1)-1);
TX_location_curr(1:NumPair/4,2) = 210*(rand(NumPair/4,1));
RX_location_curr(1:NumPair/4,1) = TX_location_curr(1:NumPair/4,1);
RX_location_curr(1:NumPair/4,2) = TX_location_curr(1:NumPair/4,2)+15;

TX_location_curr(NumPair/4+1:NumPair/2,1) = 120*ceil(3*rand(NumPair/4,1)-1);
TX_location_curr(NumPair/4+1:NumPair/2,2) = 210*(rand(NumPair/4,1))+25;
RX_location_curr(NumPair/4+1:NumPair/2,1) = TX_location_curr(NumPair/4+1:NumPair/2,1);
RX_location_curr(NumPair/4+1:NumPair/2,2) = TX_location_curr(NumPair/4+1:NumPair/2,2)-15;

TX_location_curr(NumPair/2+1:NumPair*3/4,1) = 210*(rand(NumPair/4,1));
TX_location_curr(NumPair/2+1:NumPair*3/4,2) = 120*ceil(3*rand(NumPair/4,1)-1);
RX_location_curr(NumPair/2+1:NumPair*3/4,1) = TX_location_curr(NumPair/2+1:NumPair*3/4,1)+15;
RX_location_curr(NumPair/2+1:NumPair*3/4,2) = TX_location_curr(NumPair/2+1:NumPair*3/4,2);

TX_location_curr(NumPair*3/4+1:NumPair,1) = 210*(rand(NumPair/4,1))+25;
TX_location_curr(NumPair*3/4+1:NumPair,2) = 120*ceil(3*rand(NumPair/4,1)-1);
RX_location_curr(NumPair*3/4+1:NumPair,1) = TX_location_curr(NumPair*3/4+1:NumPair,1)-15;
RX_location_curr(NumPair*3/4+1:NumPair,2) = TX_location_curr(NumPair*3/4+1:NumPair,2);

VUE_NS_direct_TX = [-1*ones(NumPair/4,1);ones(NumPair/4,1);zeros(NumPair/2,1)]; %%% +1: the north. -1: the south. 0: no
VUE_NS_direct_RX = [-1*ones(NumPair/4,1);ones(NumPair/4,1);zeros(NumPair/2,1)]; %%% +1: the north. -1: the south. 0: no
VUE_EW_direct_TX = [zeros(NumPair/2,1);-1*ones(NumPair/4,1);ones(NumPair/4,1);]; %%% +1: the east.  -1: the south. 0: no
VUE_EW_direct_RX = [zeros(NumPair/2,1);-1*ones(NumPair/4,1);ones(NumPair/4,1);]; %%% +1: the east.  -1: the south. 0: no

% plot(TX_location_curr(:,1),TX_location_curr(:,2),'bo',RX_location_curr(:,1),RX_location_curr(:,2),'rx')
% ax = gca;
% ax.XLim = [0,240];
% ax.XTick = [0:30:240];
% ax.YLim = [0,240];
% ax.YTick = [0:30:240];
%%%%%%%%%%%%%%%%%%%%%%%
avg_packet_num = Traffic_demand/avg_packet_size; %packet/s

for t = 1: simutime
    t
    Arr_pack = poissrnd(avg_packet_num*coherence_time,[NumPair,1]);
    for a = 1:NumPair
      Arrival(a,t) = sum(exprnd(avg_packet_size,[1,Arr_pack(a)]),2)/BW/coherence_time;
    end
    
   % Arrival(:,t) = exprnd(Traffic_demand/BW,[NumPair,1]);
%     Arrival(:,t) = poissrnd(Traffic_demand/BW,[NumPair,1]);
    
    %%%%%%%%%%%%%%%%%%%%%%% Spectral clustering
    if mod(t-1,T0) == 0
        mean_location = (TX_location_curr + RX_location_curr)/2;
        
        S_matrix = zeros(NumPair,NumPair);
        for a = 1 : NumPair
            for b = 1 : NumPair
                if  (norm(mean_location(a,:)-mean_location(b,:))) < delta
                    S_matrix(a,b) = exp( -( norm(mean_location(a,:)-mean_location(b,:)) )^2/(zeta^2) );
                end
            end
        end
        D_matrix = diag(sum(S_matrix,2));
        [eigV,eigD] = eig( eye(NumPair)- (D_matrix^(-0.5)) * S_matrix * (D_matrix^(-0.5)) );
        [SortValue, SortIndex]=sort( abs(transpose(diag(eigD))) , 'ascend');
        eigV = eigV(:,SortIndex);
%         D_1 = SortValue.';
%         D_2 = D_1;
%         D_1(end)=[];
%         D_2(1)=[];
%         diff_D = D_2 - D_1;
        
%         if (sum(diff_D)<0.01)
%             NumGroup = NumPair;
%         else
%             [~,NumGroup] = max(diff_D);
%         end

        U_matrix = eigV(:,[1:1:NumGroup]);
        for k = 1:NumPair
            U_matrix(k,:)=U_matrix(k,:)/norm(U_matrix(k,:));
        end
        cluster_ind = kmeans(U_matrix,NumGroup);
        
        
%         for a = 1:NumGroup
%             dummy = find(cluster_ind==a);
%             figure,plot(mean_location(dummy,1),mean_location(dummy,2),'ob');
%             name=[num2str(a) 'th group'];
%             legend(name),
%             ax = gca;
%             ax.XLim = [0,240];
%             ax.XTick = [0:30:240];
%             ax.YLim = [0,240];
%             ax.YTick = [0:30:240];
%         end
        
    end
    % %%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    %
    %
    availabel_RB_cen = zeros(NumPair, NumRB);
    for a = 1: NumGroup
        dummy = find(cluster_ind==a);
        dummy = dummy.';
        dummy2 = dummy((ceil(size(dummy,2)*rand([1,NumRB]))));
        availabel_RB_cen(dummy2+([0:1:NumRB-1]*NumPair)) =ones(1,NumRB);
    end
    availabel_RB_loc = availabel_RB_cen;
    availabel_RB_base = availabel_RB_cen;
    %
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%Path loss
    TX = repmat(TX_location_curr,[1,1,NumPair]);
    RX=[];
    RX(1,1,:) = reshape(RX_location_curr(:,1) ,[1,1,NumPair] );
    RX(1,2,:) = reshape(RX_location_curr(:,2) ,[1,1,NumPair] );
    RX = repmat(RX,[NumPair,1]);
    DX = abs(TX-RX);
    Path_loss = zeros(NumPair,NumPair);
    CHANNEL = zeros(NumPair,NumRB,NumPair);
    channel = zeros(NumPair,NumRB);
    for a=1:NumPair
        for b= 1:NumPair
            if (DX(a,1,b)==0 || DX(a,2,b)==0) && ( (TX_location_curr(a,1)==0 && RX_location_curr(b,1)==0 )|| (TX_location_curr(a,1)==120 && RX_location_curr(b,1)==120) || (TX_location_curr(a,1)==240 && RX_location_curr(b,1)==240) || (TX_location_curr(a,2)==0 && RX_location_curr(b,2)==0 )|| (TX_location_curr(a,2)==120 && RX_location_curr(b,2)==120) || (TX_location_curr(a,2)==240 && RX_location_curr(b,2)==240)  )
                Path_loss(a,b) = 10^(-68.5/10)*(max(DX(a,1,b)+DX(a,2,b),1))^(-1.61);
            elseif (min(DX(a,1,b),DX(a,2,b))<= triangle) && (min(DX(a,1,b),DX(a,2,b))> 0) && (( (TX_location_curr(a,1)==0 || TX_location_curr(a,1)==120 || TX_location_curr(a,1)==240) && (RX_location_curr(b,2)==0 || RX_location_curr(b,2)==120 || RX_location_curr(b,2)==240) ) || ( (TX_location_curr(a,2)==0 || TX_location_curr(a,2)==120 || TX_location_curr(a,2)==240) && (RX_location_curr(b,1)==0 || RX_location_curr(b,1)==120 || RX_location_curr(b,1)==240) )  )
                Path_loss(a,b) = 10^(-68.5/10)*(max(DX(a,1,b)+DX(a,2,b),1))^(-1.61);
            elseif (min(DX(a,1,b),DX(a,2,b))> triangle)  && (( (TX_location_curr(a,1)==0 || TX_location_curr(a,1)==120 || TX_location_curr(a,1)==240) && (RX_location_curr(b,2)==0 || RX_location_curr(b,2)==120 || RX_location_curr(b,2)==240) ) || ( (TX_location_curr(a,2)==0 || TX_location_curr(a,2)==120 || TX_location_curr(a,2)==240) && (RX_location_curr(b,1)==0 || RX_location_curr(b,1)==120 || RX_location_curr(b,1)==240) )  )
                Path_loss(a,b) = 10^(-54.5/10)*(DX(a,1,b)*DX(a,2,b))^(-1.61);
            end
            CHANNEL(a,:,b) = Path_loss(a,b)*exprnd(1,[1,NumRB]);
        end
        channel(a,:)= CHANNEL(a,:,a);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    %
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%Central Power allocation
    %%%%%%%%
%     dummy_queue = NetQeue_A + (2*NetQeue_B+1).*(Queue_cen(:,t)+Arrival(:,t)) + 2.*((Queue_cen(:,t)+Arrival(:,t)).^3);
%     power_cen = zeros(NumPair,NumRB);
%     for b = 1 : NumPair
%         if dummy_queue(b,1)>0
%             P_sum = P_max;
%             user_power_on_RB = zeros(1,NumRB);
%             
%             % availabel_RB %% (NumPair x NumRB)  matrix; the element of available fog is 1. Otherwise 0.
%             
%             for m = 1 : NumRB
%                 
%                 if  (availabel_RB_cen(b,m) == 1) &&  (dummy_queue(b,1)*channel(b,m)/log(2)/noise > V ) &&  (dummy_queue(b,1)*channel(b,m)/log(2)/(noise+P_sum*channel(b,m)) < V )
%                     
%                     %%%        power allocation
%                     P_do = 0;
%                     P_up = P_sum;
%                     user_power_on_RB(1,m) = P_sum/2;
%                     test_level =  dummy_queue(b,1)*channel(b,m)/log(2)/(noise+channel(b,m)*user_power_on_RB(1,m));
%                     
%                     while ( (abs(V-test_level) >= 0.002) && ((P_up-P_do)>0.001))
%                         if V < test_level
%                             P_do = user_power_on_RB(m);
%                             user_power_on_RB(m) = (P_up+P_do)/2;
%                         elseif V > test_level
%                             P_up = user_power_on_RB(m);
%                             user_power_on_RB(m) = (P_up+P_do)/2;
%                         end
%                         test_level =  dummy_queue(b,1)*channel(b,m)/log(2)/(noise+channel(b,m)*user_power_on_RB(m));
%                     end
%                 elseif  (availabel_RB_cen(b,m) == 1) &&  (dummy_queue(b,1)*channel(b,m)/log(2)/(noise+P_sum*channel(b,m)) == V )
%                     user_power_on_RB(1,m) = P_sum;
%                 elseif  (availabel_RB_cen(b,m) == 1) &&  (dummy_queue(b,1)*channel(b,m)/log(2)/(noise+P_sum*channel(b,m)) > V )
%                     user_power_on_RB(1,m) = P_sum+1;
%                 end
%                 
%             end
%             
%             %equal to P_max
%             if sum(user_power_on_RB) > P_sum
%                 
%                 noise_level = availabel_RB_cen(b,:).*channel(b,:).*dummy_queue(b,1)/log(2)/noise-V;
%                 [des_level, des_index] = sort (noise_level,'descend');
%                 
%                 cri =  availabel_RB_cen(b,des_index(1))*dummy_queue(b,1)*channel(b,des_index(1))/log(2)/(noise+P_sum*channel(b,des_index(1))) -V;
%                 des_index(find(des_level < cri))=[];
%                 
%                 if noise_level(des_index(1)) == noise_level(des_index(end))
%                     tx_P = zeros(1,NumRB);
%                     reduce_index = 0;
%                 else
%                     reduce_index = 1;
%                 end
%                 
%                 %%%Bisection
%                 
%                 while (reduce_index)
%                     cri = noise_level(des_index(end));
%                     last_index = size(des_index,2);
%                     for n = size(des_index,2)-1:-1:2
%                         if noise_level(des_index(n))== cri
%                             last_index = n;
%                         end
%                     end
%                     tx_P = zeros(1,NumRB);
%                     m=1;
%                     while (m <= last_index-1)
%                         if (dummy_queue(b,1)*channel(b,des_index(m))/log(2)/(noise+P_sum*channel(b,des_index(m)))  -V) >= cri
%                             tx_P(des_index(m))= P_sum;
%                             m = last_index;
%                         else
%                             P_up = P_sum;
%                             P_do = 0;
%                             P_test = (P_up+P_do)/2;
%                             test_level = dummy_queue(b,1)*channel(b,des_index(m))/log(2)/(noise+P_test*channel(b,des_index(m)))-V;
%                             
%                             while (abs(test_level - cri) >= 0.001)
%                                 if test_level > cri
%                                     P_do = P_test;
%                                     P_test =  (P_up+P_do)/2;
%                                 elseif test_level < cri
%                                     P_up = P_test;
%                                     P_test = (P_up+P_do)/2;
%                                 end
%                                 test_level =  dummy_queue(b,1)*channel(b,des_index(m))/log(2)/(noise+P_test*channel(b,des_index(m))) -V;
%                             end
%                             tx_P(des_index(m))= P_test;
%                             m=m+1;
%                         end
%                     end
%                     if sum(tx_P) >= P_sum
%                         des_index(last_index:end) = [];
%                     else
%                         reduce_index = 0;
%                     end
%                 end
%                 
%                 %%%        power allocation
%                 base_P = tx_P;
%                 P_last = P_sum - sum(tx_P);
%                 tx_P(des_index) = tx_P(des_index) + P_last/size(des_index,2);
%                 test_level = dummy_queue(b,1).*channel(b,des_index)/log(2)./(noise+tx_P(des_index).*channel(b,des_index))  -V;
%                 
%                 max_value = max(test_level);
%                 min_value = min(test_level);
%                 
%                 while ( (max_value-min_value >= 0.002) && (P_last>0.001))
%                     
%                     reallocate_index = des_index;
%                     reallocate_index(find(test_level>=(max_value-0.002)))=[];
%                     
%                     for m = reallocate_index
%                         P_up = tx_P(m);
%                         P_do = base_P(m);
%                         P_test = (P_up+P_do)/2;
%                         dummy_level =   dummy_queue(b,1)*channel(b,m)/log(2)/(noise+P_test*channel(b,m)) -V;
%                         
%                         while ((max_value - dummy_level)>= 0.002 || (dummy_level - max_value)>0)
%                             if dummy_level > max_value
%                                 P_do = P_test;
%                                 P_test = (P_up+P_do)/2;
%                             elseif dummy_level < max_value
%                                 P_up = P_test;
%                                 P_test = (P_up+P_do)/2;
%                             end
%                             dummy_level =  dummy_queue(b,1)*channel(b,m)/log(2)/(noise+P_test*channel(b,m)) -V;
%                         end
%                         tx_P(m)= P_test;
%                     end
%                     base_P = tx_P;
%                     P_last = P_sum - sum(tx_P);
%                     tx_P(des_index) = tx_P(des_index) + P_last/size(des_index,2);
%                     test_level =  dummy_queue(b,1).*channel(b,des_index)/log(2)./(noise+tx_P(des_index).*channel(b,des_index) ) -V;
%                     
%                     max_value = max(test_level);
%                     min_value = min(test_level);
%                 end
%                 user_power_on_RB = tx_P;
%             end
%             power_cen(b,:) = user_power_on_RB;
%         end
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    
%     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%Local power allocation
    %%%%%%%%
    dummy_queue = (LocQeue_A + (2*LocQeue_B+1).*(Queue_loc(:,t)+Arrival(:,t)) + 2.*((Queue_loc(:,t)+Arrival(:,t)).^3)).*...
        (sign(sign(1 + xi(:,t) .* (Queue_loc(:,t)+Arrival(:,t) - mu(:,t)) ./ sigma(:,t))-1)+1)+Queue_loc(:,t)+Arrival(:,t);
    
    dummy_queue(find(isnan(dummy_queue)==1)) = Queue_loc(find(isnan(dummy_queue)==1),t)+Arrival(find(isnan(dummy_queue)==1),t); 
    power_loc = zeros(NumPair,NumRB);
    
    %    active_UE = find((1 + xi(:,t) .* (Queue_loc(:,t)+Arrival(:,t) - mu(:,t)) ./ sigma(:,t))>= 0).';
    for b = 1 : NumPair
        if dummy_queue(b,1)>0
            P_sum = P_max;
            user_power_on_RB = zeros(1,NumRB);
            
            % availabel_RB %% (NumPair x NumRB)  matrix; the element of available fog is 1. Otherwise 0.
            
            for m = 1 : NumRB
                
                if  (availabel_RB_loc(b,m) == 1) &&  (dummy_queue(b,1)*channel(b,m)/log(2)/noise > V ) &&  (dummy_queue(b,1)*channel(b,m)/log(2)/(noise+P_sum*channel(b,m)) < V )
                    
                    %%%        power allocation
                    P_do = 0;
                    P_up = P_sum;
                    user_power_on_RB(1,m) = P_sum/2;
                    test_level =  dummy_queue(b,1)*channel(b,m)/log(2)/(noise+channel(b,m)*user_power_on_RB(1,m));
                    
                    while ( (abs(V-test_level) >= 0.002) && ((P_up-P_do)>0.001))
                        if V < test_level
                            P_do = user_power_on_RB(m);
                            user_power_on_RB(m) = (P_up+P_do)/2;
                        elseif V > test_level
                            P_up = user_power_on_RB(m);
                            user_power_on_RB(m) = (P_up+P_do)/2;
                        end
                        test_level =  dummy_queue(b,1)*channel(b,m)/log(2)/(noise+channel(b,m)*user_power_on_RB(m));
                    end
                elseif  (availabel_RB_loc(b,m) == 1) &&  (dummy_queue(b,1)*channel(b,m)/log(2)/(noise+P_sum*channel(b,m)) == V )
                    user_power_on_RB(1,m) = P_sum;
                elseif  (availabel_RB_loc(b,m) == 1) &&  (dummy_queue(b,1)*channel(b,m)/log(2)/(noise+P_sum*channel(b,m)) > V )
                    user_power_on_RB(1,m) = P_sum+1;
                end
                
            end
            
            %equal to P_max
            if sum(user_power_on_RB) > P_sum
                
                noise_level = availabel_RB_loc(b,:).*channel(b,:).*dummy_queue(b,1)/log(2)/noise-V;
                [des_level, des_index] = sort (noise_level,'descend');
                
                cri =  availabel_RB_loc(b,des_index(1))*dummy_queue(b,1)*channel(b,des_index(1))/log(2)/(noise+P_sum*channel(b,des_index(1))) -V;
                des_index(find(des_level < cri))=[];
                
                if noise_level(des_index(1)) == noise_level(des_index(end))
                    tx_P = zeros(1,NumRB);
                    reduce_index = 0;
                else
                    reduce_index = 1;
                end
                
                %%%Bisection
                
                while (reduce_index)
                    cri = noise_level(des_index(end));
                    last_index = size(des_index,2);
                    for n = size(des_index,2)-1:-1:2
                        if noise_level(des_index(n))== cri
                            last_index = n;
                        end
                    end
                    tx_P = zeros(1,NumRB);
                    m=1;
                    while (m <= last_index-1)
                        if (dummy_queue(b,1)*channel(b,des_index(m))/log(2)/(noise+P_sum*channel(b,des_index(m)))  -V) >= cri
                            tx_P(des_index(m))= P_sum;
                            m = last_index;
                        else
                            P_up = P_sum;
                            P_do = 0;
                            P_test = (P_up+P_do)/2;
                            test_level = dummy_queue(b,1)*channel(b,des_index(m))/log(2)/(noise+P_test*channel(b,des_index(m)))-V;
                            
                            while (abs(test_level - cri) >= 0.001)
                                if test_level > cri
                                    P_do = P_test;
                                    P_test =  (P_up+P_do)/2;
                                elseif test_level < cri
                                    P_up = P_test;
                                    P_test = (P_up+P_do)/2;
                                end
                                test_level =  dummy_queue(b,1)*channel(b,des_index(m))/log(2)/(noise+P_test*channel(b,des_index(m))) -V;
                            end
                            tx_P(des_index(m))= P_test;
                            m=m+1;
                        end
                    end
                    if sum(tx_P) >= P_sum
                        des_index(last_index:end) = [];
                    else
                        reduce_index = 0;
                    end
                end
                
                %%%        power allocation
                base_P = tx_P;
                P_last = P_sum - sum(tx_P);
                tx_P(des_index) = tx_P(des_index) + P_last/size(des_index,2);
                test_level = dummy_queue(b,1).*channel(b,des_index)/log(2)./(noise+tx_P(des_index).*channel(b,des_index))  -V;
                
                max_value = max(test_level);
                min_value = min(test_level);
                
                while ( (max_value-min_value >= 0.002) && (P_last>0.001))
                    
                    reallocate_index = des_index;
                    reallocate_index(find(test_level>=(max_value-0.002)))=[];
                    
                    for m = reallocate_index
                        P_up = tx_P(m);
                        P_do = base_P(m);
                        P_test = (P_up+P_do)/2;
                        dummy_level =   dummy_queue(b,1)*channel(b,m)/log(2)/(noise+P_test*channel(b,m)) -V;
                        
                        while ((max_value - dummy_level)>= 0.002 || (dummy_level - max_value)>0)
                            if dummy_level > max_value
                                P_do = P_test;
                                P_test = (P_up+P_do)/2;
                            elseif dummy_level < max_value
                                P_up = P_test;
                                P_test = (P_up+P_do)/2;
                            end
                            dummy_level =  dummy_queue(b,1)*channel(b,m)/log(2)/(noise+P_test*channel(b,m)) -V;
                        end
                        tx_P(m)= P_test;
                    end
                    base_P = tx_P;
                    P_last = P_sum - sum(tx_P);
                    tx_P(des_index) = tx_P(des_index) + P_last/size(des_index,2);
                    test_level =  dummy_queue(b,1).*channel(b,des_index)/log(2)./(noise+tx_P(des_index).*channel(b,des_index) ) -V;
                    
                    max_value = max(test_level);
                    min_value = min(test_level);
                end
                user_power_on_RB = tx_P;
            end
            power_loc(b,:) = user_power_on_RB;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
    
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%Baseline power allocation
%     
%     %%%%%%%%
% %     dummy_queue = BaseQeue_A + 2*Queue_base(:,t) + 2*Arrival(:,t);
%     dummy_queue = BaseQeue_A + Queue_base(:,t) + Arrival(:,t); %%% EuCNC
%     power_base = zeros(NumPair,NumRB);
%     for b = 1 : NumPair
%         if dummy_queue(b,1)>0
%             P_sum = P_max;
%             user_power_on_RB = zeros(1,NumRB);
%             
%             % availabel_RB %% (NumPair x NumRB)  matrix; the element of available fog is 1. Otherwise 0.
%             
%             for m = 1 : NumRB
%                 
%                 if  (availabel_RB_base(b,m) == 1) &&  (dummy_queue(b,1)*channel(b,m)/log(2)/noise > V ) &&  (dummy_queue(b,1)*channel(b,m)/log(2)/(noise+P_sum*channel(b,m)) < V )
%                     
%                     %%%        power allocation
%                     P_do = 0;
%                     P_up = P_sum;
%                     user_power_on_RB(1,m) = P_sum/2;
%                     test_level =  dummy_queue(b,1)*channel(b,m)/log(2)/(noise+channel(b,m)*user_power_on_RB(1,m));
%                     
%                     while ( (abs(V-test_level) >= 0.002) && ((P_up-P_do)>0.001))
%                         if V < test_level
%                             P_do = user_power_on_RB(m);
%                             user_power_on_RB(m) = (P_up+P_do)/2;
%                         elseif V > test_level
%                             P_up = user_power_on_RB(m);
%                             user_power_on_RB(m) = (P_up+P_do)/2;
%                         end
%                         test_level =  dummy_queue(b,1)*channel(b,m)/log(2)/(noise+channel(b,m)*user_power_on_RB(m));
%                     end
%                 elseif  (availabel_RB_base(b,m) == 1) &&  (dummy_queue(b,1)*channel(b,m)/log(2)/(noise+P_sum*channel(b,m)) == V )
%                     user_power_on_RB(1,m) = P_sum;
%                 elseif  (availabel_RB_base(b,m) == 1) &&  (dummy_queue(b,1)*channel(b,m)/log(2)/(noise+P_sum*channel(b,m)) > V )
%                     user_power_on_RB(1,m) = P_sum+1;
%                 end
%                 
%             end
%             
%             %equal to P_max
%             if sum(user_power_on_RB) > P_sum
%                 
%                 noise_level = availabel_RB_base(b,:).*channel(b,:).*dummy_queue(b,1)/log(2)/noise-V;
%                 [des_level, des_index] = sort (noise_level,'descend');
%                 
%                 cri =  availabel_RB_base(b,des_index(1))*dummy_queue(b,1)*channel(b,des_index(1))/log(2)/(noise+P_sum*channel(b,des_index(1))) -V;
%                 des_index(find(des_level < cri))=[];
%                 
%                 if noise_level(des_index(1)) == noise_level(des_index(end))
%                     tx_P = zeros(1,NumRB);
%                     reduce_index = 0;
%                 else
%                     reduce_index = 1;
%                 end
%                 
%                 %%%Bisection
%                 
%                 while (reduce_index)
%                     cri = noise_level(des_index(end));
%                     last_index = size(des_index,2);
%                     for n = size(des_index,2)-1:-1:2
%                         if noise_level(des_index(n))== cri
%                             last_index = n;
%                         end
%                     end
%                     tx_P = zeros(1,NumRB);
%                     m=1;
%                     while (m <= last_index-1)
%                         if (dummy_queue(b,1)*channel(b,des_index(m))/log(2)/(noise+P_sum*channel(b,des_index(m)))  -V) >= cri
%                             tx_P(des_index(m))= P_sum;
%                             m = last_index;
%                         else
%                             P_up = P_sum;
%                             P_do = 0;
%                             P_test = (P_up+P_do)/2;
%                             test_level = dummy_queue(b,1)*channel(b,des_index(m))/log(2)/(noise+P_test*channel(b,des_index(m)))-V;
%                             
%                             while (abs(test_level - cri) >= 0.001)
%                                 if test_level > cri
%                                     P_do = P_test;
%                                     P_test =  (P_up+P_do)/2;
%                                 elseif test_level < cri
%                                     P_up = P_test;
%                                     P_test = (P_up+P_do)/2;
%                                 end
%                                 test_level =  dummy_queue(b,1)*channel(b,des_index(m))/log(2)/(noise+P_test*channel(b,des_index(m))) -V;
%                             end
%                             tx_P(des_index(m))= P_test;
%                             m=m+1;
%                         end
%                     end
%                     if sum(tx_P) >= P_sum
%                         des_index(last_index:end) = [];
%                     else
%                         reduce_index = 0;
%                     end
%                 end
%                 
%                 %%%        power allocation
%                 base_P = tx_P;
%                 P_last = P_sum - sum(tx_P);
%                 tx_P(des_index) = tx_P(des_index) + P_last/size(des_index,2);
%                 test_level = dummy_queue(b,1).*channel(b,des_index)/log(2)./(noise+tx_P(des_index).*channel(b,des_index))  -V;
%                 
%                 max_value = max(test_level);
%                 min_value = min(test_level);
%                 
%                 while ( (max_value-min_value >= 0.002) && (P_last>0.001))
%                     
%                     reallocate_index = des_index;
%                     reallocate_index(find(test_level>=(max_value-0.002)))=[];
%                     
%                     for m = reallocate_index
%                         P_up = tx_P(m);
%                         P_do = base_P(m);
%                         P_test = (P_up+P_do)/2;
%                         dummy_level =   dummy_queue(b,1)*channel(b,m)/log(2)/(noise+P_test*channel(b,m)) -V;
%                         
%                         while ((max_value - dummy_level)>= 0.002 || (dummy_level - max_value)>0)
%                             if dummy_level > max_value
%                                 P_do = P_test;
%                                 P_test = (P_up+P_do)/2;
%                             elseif dummy_level < max_value
%                                 P_up = P_test;
%                                 P_test = (P_up+P_do)/2;
%                             end
%                             dummy_level =  dummy_queue(b,1)*channel(b,m)/log(2)/(noise+P_test*channel(b,m)) -V;
%                         end
%                         tx_P(m)= P_test;
%                     end
%                     base_P = tx_P;
%                     P_last = P_sum - sum(tx_P);
%                     tx_P(des_index) = tx_P(des_index) + P_last/size(des_index,2);
%                     test_level =  dummy_queue(b,1).*channel(b,des_index)/log(2)./(noise+tx_P(des_index).*channel(b,des_index) ) -V;
%                     
%                     max_value = max(test_level);
%                     min_value = min(test_level);
%                 end
%                 user_power_on_RB = tx_P;
%             end
%             power_base(b,:) = user_power_on_RB;
%         end
%     end
%     
    %%%%%%%%%%%%%%%%%%%%%%%%
%     
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%Calculate rate & update queue
%     SandI_cen = (repmat(power_cen,[1,1,NumPair])).*CHANNEL;
     SandI_loc = (repmat(power_loc,[1,1,NumPair])).*CHANNEL;
%    SandI_base = (repmat(power_base,[1,1,NumPair])).*CHANNEL;
    
    for k = 1:NumPair
%         DummySig = SandI_cen(k,:,k);
%         DummyInterfer = SandI_cen(:,:,k);
%         DummyInterfer(k,:)=[];
%         Rate_cen(k,t) =  sum(log2(1+DummySig./(noise+sum(DummyInterfer,1))),2);
   
        
         DummySig = SandI_loc(k,:,k);
         DummyInterfer = SandI_loc(:,:,k);
         DummyInterfer(k,:)=[];
         Rate_loc(k,t) =  sum(log2(1+DummySig./(noise+sum(DummyInterfer,1))),2);
   
        
%         DummySig = SandI_base(k,:,k);
%         DummyInterfer = SandI_base(:,:,k);
%         DummyInterfer(k,:)=[];
%         Rate_base(k,t) =  sum(log2(1+DummySig./(noise+sum(DummyInterfer,1))),2);
 
    end
    
%     Queue_cen(:,t+1) = max(Queue_cen(:,t) + Arrival(:,t) - Rate_cen(:,t) , zeros(NumPair,1) );
%     Actual_rate_cen(:,t) =  min(Queue_cen(:,t) + Arrival(:,t) , Rate_cen(:,t));
%     NetQeue_A =  max(NetQeue_A + max(Queue_cen(:,t+1)) - A_thr_cen, 0 );
%     NetQeue_B =  max(NetQeue_B + (max(Queue_cen(:,t+1)))^2 - B_thr_cen, 0 );
%     Mean_Pow(:,1) = Mean_Pow(:,1) + sum(power_cen,2);
%     Mean_NetQeue_A(1,t+1) = Mean_NetQeue_A(1,t) + NetQeue_A;
%     Mean_NetQeue_B(1,t+1) = Mean_NetQeue_B(1,t) + NetQeue_B;

    Queue_loc(:,t+1) = max(Queue_loc(:,t) + Arrival(:,t) - Rate_loc(:,t) , zeros(NumPair,1) );
    Actual_rate_loc(:,t) =  min(Queue_loc(:,t) + Arrival(:,t) , Rate_loc(:,t));
    LocQeue_A =  max(LocQeue_A + (Queue_loc(:,t+1) - A_thr_loc).*(sign(sign( (1+xi(:,t).*(Queue_loc(:,t+1)-mu(:,t))./sigma(:,t)) )+1)) , zeros(NumPair,1) );
    LocQeue_B =  max(LocQeue_B + ((Queue_loc(:,t+1)).^2 - B_thr_loc).*(sign(sign( (1+xi(:,t).*(Queue_loc(:,t+1)-mu(:,t))./sigma(:,t)) )+1)) , zeros(NumPair,1)  );
    Mean_Pow(:,2) = Mean_Pow(:,2) + sum(power_loc,2);
    Mean_LocQeue_A(:,t+1) = Mean_LocQeue_A(:,t) + LocQeue_A;
    Mean_LocQeue_B(:,t+1) = Mean_LocQeue_B(:,t) + LocQeue_B;
        
    
    for k = 1:NumPair
                
        [ycdf,xcdf] = cdfcalc(Queue_loc(k,:));
        dummy = find(ycdf>=1-1/NumPair);
        mu(k,t+1)=xcdf(dummy(1)-1);
        
        dummy = find(ycdf >= 1 - POT_upper);
        d = xcdf(dummy(1)-1);
        
        C_mean = mean((Queue_loc(k, find(Queue_loc(k,:)>d) )-d),2);
        C_var = mean((Queue_loc(k, find(Queue_loc(k,:)>d) )-d).^2,2);
        

       if isnan(C_mean) || isnan(C_var) || (C_var==C_mean^2)
           sigma(k,t+1)=0;
           xi(k,t+1) =0;
       else
           sigma(k,t+1) = ( C_var*C_mean + (C_var-2*C_mean^2)*(mu(k,t+1)-d) )/2/(C_var-C_mean^2);
           xi(k,t+1) = (C_var-2*C_mean^2)/2/(C_var-C_mean^2);
       end
    end
    
    
%     Queue_base(:,t+1) = max(Queue_base(:,t) + Arrival(:,t) - Rate_base(:,t) , zeros(NumPair,1) );
%     Actual_rate_base(:,t) =  min(Queue_base(:,t) + Arrival(:,t) , Rate_base(:,t));
% %     BaseQeue_A = max(BaseQeue_A + (sign(sign(Queue_base(:,t+1)-Q_th)-1)+1)- viola_prob, zeros(NumPair,1) );
%         BaseQeue_A = max(BaseQeue_A + Queue_base(:,t+1)- A_thr_loc, zeros(NumPair,1) ); %%%EuCNC
%     Mean_BaseQeue_A(:,t+1) = Mean_BaseQeue_A(:,t) + BaseQeue_A;
%     Mean_Pow(:,3) = Mean_Pow(:,3) + sum(power_base,2);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    %
    %
    %
    %
    %
    %
    %%%%%%%%%%%% Updating location
    instan_speed_TX = Avg_speed + (acc_speed*rand(NumPair,1)-acc_speed/2)*coherence_time; %%% meter/sec
    instan_speed_RX = Avg_speed + (acc_speed*rand(NumPair,1)-acc_speed/2)*coherence_time; %%% meter/sec
    location_increase_TX = coherence_time * instan_speed_TX;
    location_increase_RX = coherence_time * instan_speed_RX;
    
    
    TX_location_next(:,1) =  TX_location_curr(:,1) + VUE_EW_direct_TX.* location_increase_TX;
    TX_location_next(:,2) =  TX_location_curr(:,2) + VUE_NS_direct_TX.* location_increase_TX;
    RX_location_next(:,1) =  RX_location_curr(:,1) + VUE_EW_direct_RX.* location_increase_RX;
    RX_location_next(:,2) =  RX_location_curr(:,2) + VUE_NS_direct_RX.* location_increase_RX;
    %%%%%%%%%%%%%%%%%
    
    % %%%%%%%%%%%%%%%%% Updating while near intersection for TX
    dummy_ind_1 = (find( ( (TX_location_next(:,1)    ) .* (TX_location_curr(:,1)    ) ) < 0 )).';
    dummy_ind_2 = (find( ( (TX_location_next(:,1)-120) .* (TX_location_curr(:,1)-120) ) < 0 )).';
    dummy_ind_3 = (find( ( (TX_location_next(:,1)-240) .* (TX_location_curr(:,1)-240) ) < 0 )).';
    dummy_ind_4 = (find( ( (TX_location_next(:,2)    ) .* (TX_location_curr(:,2)    ) ) < 0 )).';
    dummy_ind_5 = (find( ( (TX_location_next(:,2)-120) .* (TX_location_curr(:,2)-120) ) < 0 )).';
    dummy_ind_6 = (find( ( (TX_location_next(:,2)-240) .* (TX_location_curr(:,2)-240) ) < 0 )).';
    
    for k = dummy_ind_1
        if TX_location_curr(k,2)==0
            TX_location_next(k,1) = 0;
            TX_location_next(k,2) = location_increase_TX(k) - TX_location_curr(k,1);
            VUE_EW_direct_TX(k,1) = 0;
            VUE_NS_direct_TX(k,1) = 1;
        elseif TX_location_curr(k,2)== 120
            TX_location_next(k,1) = 0;
            dummy_direction = 2*round(rand(1))-1;
            TX_location_next(k,2) = 120 + dummy_direction*(location_increase_TX(k) - TX_location_curr(k,1));
            VUE_EW_direct_TX(k,1) = 0;
            VUE_NS_direct_TX(k,1) = dummy_direction;
        elseif TX_location_curr(k,2)== 240
            TX_location_next(k,1) = 0;
            TX_location_next(k,2) = 240 - (location_increase_TX(k) - TX_location_curr(k,1));
            VUE_EW_direct_TX(k,1) = 0;
            VUE_NS_direct_TX(k,1) = -1;
        end
    end
    
    
    
    for k = dummy_ind_2
        if (TX_location_curr(k,1)< 120)
            if TX_location_curr(k,2)==0
                dummy_direction = round(rand(1));
                VUE_EW_direct_TX(k,1) = mod(VUE_EW_direct_TX(k,1)+dummy_direction,2);
                VUE_NS_direct_TX(k,1) = dummy_direction;
                reman = location_increase_TX(k) + TX_location_curr(k,1) - 120;
                TX_location_next(k,1) = 120 + VUE_EW_direct_TX(k,1)*reman;
                TX_location_next(k,2) =  VUE_NS_direct_TX(k,1)*reman;
                
            elseif TX_location_curr(k,2)== 120
                dummy_direction = ceil(3*rand(1))-2;
                if (dummy_direction==1) || (dummy_direction == -1)
                    VUE_EW_direct_TX(k,1) = 0;
                    VUE_NS_direct_TX(k,1) = dummy_direction;
                end
                reman = location_increase_TX(k) + TX_location_curr(k,1)-120;
                TX_location_next(k,1) = 120 + VUE_EW_direct_TX(k,1)*reman;
                TX_location_next(k,2) = 120 + VUE_NS_direct_TX(k,1)*reman;
                
            elseif TX_location_curr(k,2) == 240
                dummy_direction = round(rand(1))-1;
                VUE_EW_direct_TX(k,1) = mod(VUE_EW_direct_TX(k,1)+dummy_direction,2);
                VUE_NS_direct_TX(k,1) = dummy_direction;
                reman = location_increase_TX(k) + TX_location_curr(k,1)-120;
                TX_location_next(k,1) = 120 + VUE_EW_direct_TX(k,1)*reman;
                TX_location_next(k,2) =  240 + VUE_NS_direct_TX(k,1)*reman;
            end
        elseif (TX_location_curr(k,1)> 120)
            if TX_location_curr(k,2)==0
                dummy_direction = round(rand(1));
                VUE_EW_direct_TX(k,1) = VUE_EW_direct_TX(k,1)+dummy_direction;
                VUE_NS_direct_TX(k,1) = dummy_direction;
                reman = location_increase_TX(k) - (TX_location_curr(k,1)-120);
                TX_location_next(k,1) = 120 + VUE_EW_direct_TX(k,1)*reman;
                TX_location_next(k,2) =  VUE_NS_direct_TX(k,1)*reman;
                
            elseif TX_location_curr(k,2)== 120
                dummy_direction = ceil(3*rand(1))-2;
                if (dummy_direction==1) || (dummy_direction == -1)
                    VUE_EW_direct_TX(k,1) = 0;
                    VUE_NS_direct_TX(k,1) = dummy_direction;
                end
                reman = location_increase_TX(k) - TX_location_curr(k,1)+120;
                TX_location_next(k,1) = 120 + VUE_EW_direct_TX(k,1)*reman;
                TX_location_next(k,2) = 120 + VUE_NS_direct_TX(k,1)*reman;
                
            elseif TX_location_curr(k,2)== 240
                dummy_direction = round(rand(1))-1;
                VUE_EW_direct_TX(k,1) = mod(VUE_EW_direct_TX(k,1)+dummy_direction,-2);
                VUE_NS_direct_TX(k,1) = dummy_direction;
                reman = location_increase_TX(k) - TX_location_curr(k,1)+120;
                TX_location_next(k,1) = 120 + VUE_EW_direct_TX(k,1)*reman;
                TX_location_next(k,2) = 240 + VUE_NS_direct_TX(k,1)*reman;
            end
        end
    end
    
    
    
    
    for k = dummy_ind_3
        if TX_location_curr(k,2)==0
            TX_location_next(k,1) = 240;
            TX_location_next(k,2) = location_increase_TX(k) - (240-TX_location_curr(k,1));
            VUE_EW_direct_TX(k,1) = 0;
            VUE_NS_direct_TX(k,1) = 1;
        elseif TX_location_curr(k,2)== 120
            TX_location_next(k,1) = 240;
            dummy_direction = 2*round(rand(1))-1;
            TX_location_next(k,2) = 120 + dummy_direction*(location_increase_TX(k) - (240-TX_location_curr(k,1)) );
            VUE_EW_direct_TX(k,1) = 0;
            VUE_NS_direct_TX(k,1) = dummy_direction;
        elseif TX_location_curr(k,2)== 240
            TX_location_next(k,1) = 240;
            TX_location_next(k,2) = 240 - (location_increase_TX(k) - (240-TX_location_curr(k,1)));
            VUE_EW_direct_TX(k,1) = 0;
            VUE_NS_direct_TX(k,1) = -1;
        end
    end
    
    
    
    
    
    for k = dummy_ind_4
        if TX_location_curr(k,1)==0
            TX_location_next(k,2) = 0;
            TX_location_next(k,1) = location_increase_TX(k) - TX_location_curr(k,2);
            VUE_EW_direct_TX(k,1) = 1;
            VUE_NS_direct_TX(k,1) = 0;
        elseif TX_location_curr(k,1)== 120
            TX_location_next(k,2) = 0;
            dummy_direction = 2*round(rand(1))-1;
            TX_location_next(k,1) = 120 + dummy_direction*(location_increase_TX(k) - TX_location_curr(k,2));
            VUE_EW_direct_TX(k,1) = dummy_direction;
            VUE_NS_direct_TX(k,1) = 0;
        elseif TX_location_curr(k,1)== 240
            TX_location_next(k,2) = 0;
            TX_location_next(k,1) = 240 - (location_increase_TX(k) - TX_location_curr(k,2));
            VUE_EW_direct_TX(k,1) = -1;
            VUE_NS_direct_TX(k,1) = 0;
        end
    end
    
    
    for k = dummy_ind_5
        if (TX_location_curr(k,2)< 120)
            if TX_location_curr(k,1)==0
                dummy_direction = round(rand(1));
                VUE_NS_direct_TX(k,1) = mod(VUE_NS_direct_TX(k,1)+dummy_direction,2);
                VUE_EW_direct_TX(k,1) = dummy_direction;
                reman = location_increase_TX(k) + TX_location_curr(k,2) - 120;
                TX_location_next(k,1) =  VUE_EW_direct_TX(k,1)*reman;
                TX_location_next(k,2) = 120 + VUE_NS_direct_TX(k,1)*reman;
                
            elseif TX_location_curr(k,1)== 120
                dummy_direction = ceil(3*rand(1))-2;
                if (dummy_direction==1) || (dummy_direction == -1)
                    VUE_NS_direct_TX(k,1) = 0;
                    VUE_EW_direct_TX(k,1) = dummy_direction;
                end
                reman = location_increase_TX(k) + TX_location_curr(k,2)-120;
                TX_location_next(k,1) = 120 + VUE_EW_direct_TX(k,1)*reman;
                TX_location_next(k,2) = 120 + VUE_NS_direct_TX(k,1)*reman;
                
            elseif TX_location_curr(k,1)== 240
                dummy_direction = round(rand(1))-1;
                VUE_NS_direct_TX(k,1) = mod(VUE_NS_direct_TX(k,1)+dummy_direction,2);
                VUE_EW_direct_TX(k,1) = dummy_direction;
                reman = location_increase_TX(k) + TX_location_curr(k,2)-120;
                TX_location_next(k,1) = 240 + VUE_EW_direct_TX(k,1)*reman;
                TX_location_next(k,2) =  120 + VUE_NS_direct_TX(k,1)*reman;
            end
            
        elseif (TX_location_curr(k,2)> 120)
            
            if TX_location_curr(k,1)==0
                dummy_direction = round(rand(1));
                VUE_NS_direct_TX(k,1) = VUE_NS_direct_TX(k,1)+dummy_direction;
                VUE_EW_direct_TX(k,1) = dummy_direction;
                reman = location_increase_TX(k) - (TX_location_curr(k,2)-120);
                TX_location_next(k,1) =  VUE_EW_direct_TX(k,1)*reman;
                TX_location_next(k,2) = 120 + VUE_NS_direct_TX(k,1)*reman;
                
            elseif TX_location_curr(k,1)== 120
                dummy_direction = ceil(3*rand(1))-2;
                if (dummy_direction==1) || (dummy_direction == -1)
                    VUE_NS_direct_TX(k,1) = 0;
                    VUE_EW_direct_TX(k,1) = dummy_direction;
                end
                reman = location_increase_TX(k) - TX_location_curr(k,2)+120;
                TX_location_next(k,1) = 120 + VUE_EW_direct_TX(k,1)*reman;
                TX_location_next(k,2) = 120 + VUE_NS_direct_TX(k,1)*reman;
                
            elseif TX_location_curr(k,1)== 240
                dummy_direction = round(rand(1))-1;
                VUE_NS_direct_TX(k,1) = mod(VUE_NS_direct_TX(k,1)+dummy_direction,-2);
                VUE_EW_direct_TX(k,1) = dummy_direction;
                reman = location_increase_TX(k) - TX_location_curr(k,2)+120;
                TX_location_next(k,1) = 240 + VUE_EW_direct_TX(k,1)*reman;
                TX_location_next(k,2) =  120 + VUE_NS_direct_TX(k,1)*reman;
            end
        end
    end
    
    
    
    for k = dummy_ind_6
        if TX_location_curr(k,1)==0
            TX_location_next(k,2) = 240;
            TX_location_next(k,1) = location_increase_TX(k) - (240-TX_location_curr(k,2));
            VUE_EW_direct_TX(k,1) = 1;
            VUE_NS_direct_TX(k,1) = 0;
        elseif TX_location_curr(k,1)== 120
            TX_location_next(k,2) = 240;
            dummy_direction = 2*round(rand(1))-1;
            TX_location_next(k,1) = 120 + dummy_direction*(location_increase_TX(k) - (240-TX_location_curr(k,2)) );
            VUE_NS_direct_TX(k,1) = 0;
            VUE_EW_direct_TX(k,1) = dummy_direction;
        elseif TX_location_curr(k,1)== 240
            TX_location_next(k,2) = 240;
            TX_location_next(k,1) = 240 - (location_increase_TX(k) - (240-TX_location_curr(k,2)) );
            VUE_EW_direct_TX(k,1) = -1;
            VUE_NS_direct_TX(k,1) = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Updating while near intersection for RX
    dummy_ind_1 = (find( ( (RX_location_next(:,1)    ) .* (RX_location_curr(:,1)    ) ) < 0 )).';
    dummy_ind_2 = (find( ( (RX_location_next(:,1)-120) .* (RX_location_curr(:,1)-120) ) < 0 )).';
    dummy_ind_3 = (find( ( (RX_location_next(:,1)-240) .* (RX_location_curr(:,1)-240) ) < 0 )).';
    dummy_ind_4 = (find( ( (RX_location_next(:,2)    ) .* (RX_location_curr(:,2)    ) ) < 0 )).';
    dummy_ind_5 = (find( ( (RX_location_next(:,2)-120) .* (RX_location_curr(:,2)-120) ) < 0 )).';
    dummy_ind_6 = (find( ( (RX_location_next(:,2)-240) .* (RX_location_curr(:,2)-240) ) < 0 )).';
    
    
    
    for k = dummy_ind_1
        if RX_location_curr(k,2)==0
            RX_location_next(k,1) = 0;
            RX_location_next(k,2) = location_increase_RX(k) - RX_location_curr(k,1);
            VUE_EW_direct_RX(k,1) = 0;
            VUE_NS_direct_RX(k,1) = 1;
        elseif RX_location_curr(k,2)== 120
            RX_location_next(k,1) = 0;
            RX_location_next(k,2) = 120 + VUE_NS_direct_TX(k,1)*(location_increase_RX(k) - RX_location_curr(k,1));
            VUE_EW_direct_RX(k,1) = 0;
            VUE_NS_direct_RX(k,1) = VUE_NS_direct_TX(k,1);
            
        elseif RX_location_curr(k,2)== 240
            RX_location_next(k,1) = 0;
            RX_location_next(k,2) = 240 - (location_increase_RX(k) - RX_location_curr(k,1));
            VUE_EW_direct_RX(k,1) = 0;
            VUE_NS_direct_RX(k,1) = -1;
        end
    end
    
    
    
    for k = dummy_ind_2
        if (RX_location_curr(k,1)< 120)
            if RX_location_curr(k,2)==0
                VUE_EW_direct_RX(k,1) = VUE_EW_direct_TX(k,1);
                VUE_NS_direct_RX(k,1) = VUE_NS_direct_TX(k,1);
                reman = location_increase_RX(k) + RX_location_curr(k,1) - 120;
                RX_location_next(k,1) = 120 + VUE_EW_direct_RX(k,1)*reman;
                RX_location_next(k,2) =  VUE_NS_direct_RX(k,1)*reman;
                
            elseif RX_location_curr(k,2)== 120
                
                VUE_EW_direct_RX(k,1) = VUE_EW_direct_TX(k,1);
                VUE_NS_direct_RX(k,1) = VUE_NS_direct_TX(k,1);
                reman = location_increase_RX(k) + RX_location_curr(k,1)-120;
                RX_location_next(k,1) = 120 + VUE_EW_direct_RX(k,1)*reman;
                RX_location_next(k,2) = 120 + VUE_NS_direct_RX(k,1)*reman;
                
            elseif RX_location_curr(k,2)== 240
                
                VUE_EW_direct_RX(k,1) = VUE_EW_direct_TX(k,1);
                VUE_NS_direct_RX(k,1) = VUE_NS_direct_TX(k,1);
                reman = location_increase_RX(k) + RX_location_curr(k,1)-120;
                RX_location_next(k,1) = 120 + VUE_EW_direct_RX(k,1)*reman;
                RX_location_next(k,2) =  240 + VUE_NS_direct_RX(k,1)*reman;
            end
        elseif (RX_location_curr(k,1)> 120)
            if RX_location_curr(k,2)==0
                VUE_EW_direct_RX(k,1) = VUE_EW_direct_TX(k,1);
                VUE_NS_direct_RX(k,1) = VUE_NS_direct_TX(k,1);
                reman = location_increase_RX(k) - (RX_location_curr(k,1)-120);
                RX_location_next(k,1) = 120 + VUE_EW_direct_RX(k,1)*reman;
                RX_location_next(k,2) =  VUE_NS_direct_RX(k,1)*reman;
                
            elseif RX_location_curr(k,2)== 120
                VUE_EW_direct_RX(k,1) = VUE_EW_direct_TX(k,1);
                VUE_NS_direct_RX(k,1) = VUE_NS_direct_TX(k,1);
                reman = location_increase_RX(k) - RX_location_curr(k,1)+120;
                RX_location_next(k,1) = 120 + VUE_EW_direct_RX(k,1)*reman;
                RX_location_next(k,2) = 120 + VUE_NS_direct_RX(k,1)*reman;
                
            elseif RX_location_curr(k,2)== 240
                VUE_EW_direct_RX(k,1) = VUE_EW_direct_TX(k,1);
                VUE_NS_direct_RX(k,1) = VUE_NS_direct_TX(k,1);
                reman = location_increase_RX(k) - RX_location_curr(k,1)+120;
                RX_location_next(k,1) = 120 + VUE_EW_direct_RX(k,1)*reman;
                RX_location_next(k,2) =  240 + VUE_NS_direct_RX(k,1)*reman;
            end
        end
    end
    
    
    
    
    for k = dummy_ind_3
        if RX_location_curr(k,2)==0
            RX_location_next(k,1) = 240;
            RX_location_next(k,2) = location_increase_RX(k) - (240-RX_location_curr(k,1));
            VUE_EW_direct_RX(k,1) = 0;
            VUE_NS_direct_RX(k,1) = 1;
        elseif RX_location_curr(k,2)== 120
            RX_location_next(k,1) = 240;
            VUE_EW_direct_RX(k,1) = 0;
            VUE_NS_direct_RX(k,1) = VUE_NS_direct_TX(k,1);
            RX_location_next(k,2) = 120 + VUE_NS_direct_RX(k,1)*(location_increase_RX(k) - (240-RX_location_curr(k,1)) );
            
        elseif RX_location_curr(k,2)== 240
            RX_location_next(k,1) = 240;
            RX_location_next(k,2) = 240 - (location_increase_RX(k) - (240-RX_location_curr(k,1)));
            VUE_EW_direct_RX(k,1) = 0;
            VUE_NS_direct_RX(k,1) = -1;
        end
    end
    
    
    
    
    
    for k = dummy_ind_4
        if RX_location_curr(k,1)==0
            RX_location_next(k,2) = 0;
            RX_location_next(k,1) = location_increase_RX(k) - RX_location_curr(k,2);
            VUE_EW_direct_RX(k,1) = 1;
            VUE_NS_direct_RX(k,1) = 0;
        elseif RX_location_curr(k,1)== 120
            VUE_EW_direct_RX(k,1) = VUE_EW_direct_TX(k,1);
            VUE_NS_direct_RX(k,1) = 0;
            RX_location_next(k,2) = 0;
            RX_location_next(k,1) = 120 + VUE_EW_direct_RX(k,1) *(location_increase_RX(k) - RX_location_curr(k,2));
            
        elseif RX_location_curr(k,1)== 240
            RX_location_next(k,2) = 0;
            RX_location_next(k,1) = 240 - (location_increase_RX(k) - RX_location_curr(k,2));
            VUE_EW_direct_RX(k,1) = -1;
            VUE_NS_direct_RX(k,1) = 0;
        end
    end
    
    
    for k = dummy_ind_5
        if RX_location_curr(k,2)< 120
            
            if RX_location_curr(k,1)==0
                VUE_NS_direct_RX(k,1) = VUE_NS_direct_TX(k,1);
                VUE_EW_direct_RX(k,1) = VUE_EW_direct_TX(k,1);
                reman = location_increase_RX(k) + RX_location_curr(k,2) - 120;
                RX_location_next(k,1) =  VUE_EW_direct_RX(k,1)*reman;
                RX_location_next(k,2) = 120 + VUE_NS_direct_RX(k,1)*reman;
            elseif RX_location_curr(k,1)== 120
                VUE_NS_direct_RX(k,1) = VUE_NS_direct_TX(k,1);
                VUE_EW_direct_RX(k,1) = VUE_EW_direct_TX(k,1);
                reman = location_increase_RX(k) + RX_location_curr(k,2)-120;
                RX_location_next(k,1) = 120 + VUE_EW_direct_RX(k,1)*reman;
                RX_location_next(k,2) = 120 + VUE_NS_direct_RX(k,1)*reman;
                
            elseif RX_location_curr(k,1)== 240
                VUE_NS_direct_RX(k,1) = VUE_NS_direct_TX(k,1);
                VUE_EW_direct_RX(k,1) = VUE_EW_direct_TX(k,1);
                reman = location_increase_RX(k) + RX_location_curr(k,2)-120;
                RX_location_next(k,1) = 240 + VUE_EW_direct_RX(k,1)*reman;
                RX_location_next(k,2) =  120 + VUE_NS_direct_RX(k,1)*reman;
            end
            
        elseif RX_location_curr(k,2)> 120
            if RX_location_curr(k,1)==0
                VUE_NS_direct_RX(k,1) = VUE_NS_direct_TX(k,1);
                VUE_EW_direct_RX(k,1) = VUE_EW_direct_TX(k,1);
                reman = location_increase_RX(k) - (RX_location_curr(k,2)-120);
                RX_location_next(k,1) =  VUE_EW_direct_RX(k,1)*reman;
                RX_location_next(k,2) = 120 + VUE_NS_direct_RX(k,1)*reman;
                
            elseif RX_location_curr(k,1)== 120
                VUE_NS_direct_RX(k,1) = VUE_NS_direct_TX(k,1);
                VUE_EW_direct_RX(k,1) = VUE_EW_direct_TX(k,1);
                reman = location_increase_RX(k) - RX_location_curr(k,2)+120;
                RX_location_next(k,1) = 120 + VUE_EW_direct_RX(k,1)*reman;
                RX_location_next(k,2) = 120 + VUE_NS_direct_RX(k,1)*reman;
                
            elseif RX_location_curr(k,1)== 240
                VUE_NS_direct_RX(k,1) = VUE_NS_direct_TX(k,1) ;
                VUE_EW_direct_RX(k,1) =  VUE_EW_direct_TX(k,1);
                reman = location_increase_RX(k) - RX_location_curr(k,2)+120;
                RX_location_next(k,1) = 240 + VUE_EW_direct_RX(k,1)*reman;
                RX_location_next(k,2) =  120 + VUE_NS_direct_RX(k,1)*reman;
            end
        end
    end
    
    
    
    for k = dummy_ind_6
        if RX_location_curr(k,1)==0
            RX_location_next(k,2) = 240;
            RX_location_next(k,1) = location_increase_RX(k) - (240-RX_location_curr(k,2));
            VUE_EW_direct_RX(k,1) = 1;
            VUE_NS_direct_RX(k,1) = 0;
        elseif RX_location_curr(k,1)== 120
            RX_location_next(k,2) = 240;
            VUE_NS_direct_RX(k,1) = VUE_NS_direct_TX(k,1);
            VUE_EW_direct_RX(k,1) = VUE_EW_direct_TX(k,1);
            RX_location_next(k,1) = 120 + VUE_EW_direct_RX(k,1)*(location_increase_RX(k) - (240-RX_location_curr(k,2)) );
            
        elseif RX_location_curr(k,1)== 240
            RX_location_next(k,2) = 240;
            RX_location_next(k,1) = 240 - (location_increase_RX(k) - (240-RX_location_curr(k,2)));
            VUE_EW_direct_RX(k,1) = -1;
            VUE_NS_direct_RX(k,1) = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    %
    %
    %
    TX_location_curr = TX_location_next;
    RX_location_curr = RX_location_next;
%     plot(TX_location_curr(:,1),TX_location_curr(:,2),'bo',RX_location_curr(:,1),RX_location_curr(:,2),'rx')
%     ax = gca;
%     ax.XLim = [0,240];
%     ax.XTick = [0:30:240];
%     ax.YLim = [0,240];
%     ax.YTick = [0:30:240];
%     tic,pause(0.0001),toc
end %%% t-loop for simulation timeline

% figure,plot([1:simutime],Mean_NetQeue_A(2:simutime+1)./[1:simutime]);
% legend('Mean NetQeue A')
% figure,plot([1:simutime],Mean_NetQeue_B(2:simutime+1)./[1:simutime]);
% legend('Mean NetQeue B')
% figure,plot([1:simutime],Mean_LocQeue_A(:,2:simutime+1)./(repmat([1:simutime],[NumPair,1])));
% legend('Mean LocQeue A')
% figure,plot([1:simutime],Mean_LocQeue_B(:,2:simutime+1)./(repmat([1:simutime],[NumPair,1])));
% legend('Mean LocQeue B')
figure,plot([1:simutime],Mean_BaseQeue_A(:,2:simutime+1)./(repmat([1:simutime],[NumPair,1])));
legend('Mean BaseQeue A')
Mean_Pow = Mean_Pow/simutime;
 name=['V2V_letter_N=' num2str(NumPair) 'V=' num2str(V)];
 save(name,'Mean_Pow','Rate_loc','Rate_base','Actual_rate_loc','Actual_rate_base','Arrival','Queue_loc','Queue_base','mu','sigma','xi')
%name=['V2V_letter_N=' num2str(NumPair) 'V=' num2str(V)];
%save(name,'Mean_Pow','Rate_base','Actual_rate_base','Arrival','Queue_base')
