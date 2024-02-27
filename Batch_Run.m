% Batch_Run
Clear_data_Keep_debug;
disp('Batch run started.');

% Set_parameters;
% VUE = 20:20:200;
% 
% for batchIndex = 1:length(VUE)
%     disp('####################################');
%     disp(['############## ' num2str( VUE(batchIndex), '%04d') ' ################']);
%     disp('####################################');
%     clearvars -except VUE batchIndex
%     
%     Set_parameters;
%     NetworkParams.TotalVUEPairs = VUE(batchIndex);
%     
%     Mainfile;
%     
% end

% start_topology = 1;
for test_topology = fliplr(0:100)
    if exist(['Completed_' num2str(test_topology,'%03d') '.mat'],'file')
        break;
    end
end
start_topology = test_topology + 1;
% start_clock = clock;
for TOPOLOGY = start_topology:100
    disp(['Topology: ' num2str(TOPOLOGY)]);
    RANDOM_SEED = TOPOLOGY;
    %%
    Set_Parameters;
    Mainfile;
    disp(['Topology ' num2str(TOPOLOGY) ' is done.']);
    pause(3);
    clearvars -except TOPOLOGY 
    %% delete and create
    if exist(['Completed_' num2str(TOPOLOGY-1,'%03d') '.mat'],'file')
        delete(['Completed_' num2str(TOPOLOGY-1,'%03d') '.mat']);
        pause(3);
    end
    save(['Completed_' num2str(TOPOLOGY,'%03d') '.mat'], 'TOPOLOGY' );
end

disp('Batch run finished.');