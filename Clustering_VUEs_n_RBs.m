function [...
    Resource_allocation_table ...
    ] = Clustering_VUEs_n_RBs( ...
    queues_now, ...
    TX_location_shifted, ...
    RX_location_shifted ...
    )

global total_VUE_pairs total_resource_blocks 
global total_clusters neighborhood similarity_std

%% Spectral clustering



%% VUEs
mean_location = (TX_location_shifted + RX_location_shifted)/2;

similarity_matrix = zeros( total_VUE_pairs, total_VUE_pairs);
for pair_a = 1 : total_VUE_pairs
    for pair_b = 1 : total_VUE_pairs
        if  (norm(mean_location(pair_a,:)-mean_location(pair_b,:))) ...
                < neighborhood
            similarity_matrix(pair_a,pair_b) = ...
                exp( -( norm(mean_location(pair_a,:)-...
                mean_location(pair_b,:)) )^2/...
                (2*similarity_std^2) );
        end
    end
end
diagonal_matrix = diag(sum(similarity_matrix,2));
% EVD
laplacian_matrix = eye(total_VUE_pairs)- (diagonal_matrix^(-0.5)) ...
    * similarity_matrix * (diagonal_matrix^(-0.5));
% Matlab has precison errors. In which the symmetricity of laplacian cannot
% be obtained and thus, eigenvalues become complex. Following line is to
% fix that issue.
laplacian_matrix = .5*( laplacian_matrix + laplacian_matrix' );

[eigVals,eigD] = eig( laplacian_matrix );
[eigSortValue, eigSortIndex]=sort( abs(transpose(diag(eigD))) , 'ascend');
eigVals = eigVals(:,eigSortIndex);
%         D_1 = SortValue.';
%         D_2 = D_1;
%         D_1(end)=[];
%         D_2(1)=[];
%         diff_D = D_2 - D_1;

%         if (sum(diff_D)<0.01)
%             NumGroup = total_VUE_pairs;
%         else
%             [~,NumGroup] = max(diff_D);
%         end

projected_locations_matrix = eigVals(:,[1:1:total_clusters]);
for pair = 1:total_VUE_pairs
    projected_locations_matrix(pair,:) = ...
        projected_locations_matrix(pair,:)/...
        norm(projected_locations_matrix(pair,:));
end
cluster_index = kmeans(projected_locations_matrix,total_clusters);


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


% %%%%%%%%%%%%%%%%%%%%%%%%
%


%% RBs

Resource_allocation_table = zeros(total_VUE_pairs, total_resource_blocks);
offset_queue = 10;
for this_cluster_index = 1: total_clusters
    this_cluster = find(cluster_index==this_cluster_index);
    this_cluster = this_cluster.';
    member_queues = queues_now(this_cluster) + offset_queue;
    member_RB_count = ...
        round( total_resource_blocks*member_queues/sum(member_queues) );
    % Need to test following for posibilitites on looping forever
    if ( sum(member_RB_count)>total_resource_blocks )
        member_RB_count = ...
            floor( total_resource_blocks*member_RB_count/sum(member_RB_count) );
    end
    if ( sum(member_RB_count)<total_resource_blocks )
        members_favored = randperm( length(this_cluster), ...
            total_resource_blocks - sum(member_RB_count));
        member_RB_count(members_favored) = ...
            member_RB_count(members_favored) + 1;
    end
    members_per_RB_ordered = [];
    for member_index = 1:length(this_cluster)
        this_member = this_cluster(member_index);
        this_RBs = member_RB_count(member_index);
        members_per_RB_ordered = ...
            [members_per_RB_ordered this_member*ones(1, this_RBs)];
    end
    members_per_RB = members_per_RB_ordered(randperm(total_resource_blocks));
    % members_per_RB = this_cluster((ceil(size(this_cluster,2)*rand([1,total_resource_blocks]))));
    Resource_allocation_table(...
        members_per_RB+([0:1:total_resource_blocks-1]*total_VUE_pairs)) ...
        = ones(1,total_resource_blocks);
end


end