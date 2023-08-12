
clear

% Perform community detection
addpath(genpath('/home/mgell/Work/Impulsivity_networks/code/'));
addpath(genpath('/home/mgell/Matlab_toolboxes/HierarchicalConsensus-master'));
addpath(genpath('/home/mgell/Matlab_toolboxes/GenLouvain-2.2.0'));

n_perm = 1000
nodes = 0       % 'extended'

% FC
W = csvread('/home/mgell/Work/Impulsivity_networks/res/all_mean_s2sFC.csv');
% without close neighbours
W([3,17,24],:) = [];
W(:,[3,17,24]) = [];

% Nodes
tab = readtable('/home/mgell/Work/Impulsivity_networks/nets/All_localmax.node', 'FileType','text');
% without close neighbours
tab([3,17,24],:) = [];


% normalise by euclidian distance
% node_coords = table2array(tab(:,1:3));
% dist_mat = coord2dist(node_coords);
% W = dist_mat./W;

% cluster
% Iterative community finetuning.
% W is the input connection matrix.
S = zeros(length(W),n_perm);
for part_n = 1:n_perm
    n  = size(W,1);             % number of nodes
    M  = 1:n;                   % initial community affiliations
    Q0 = -1; Q1 = 0;            % initialize modularity values
    while Q1-Q0>1e-5;           % while modularity increases
        Q0 = Q1;                % perform community detection
        [M, Q1] = community_louvain(W, [], M, 'negative_asym');
    end
    S(:,part_n) = M;
end
clear Q0 Q1


% Consesnsus clustering
%[Sc, Tree] = hierarchicalConsensus(S);
%drawHierarchy(Sc, Tree)

C = coclassificationMatrix(S);
%consensusPlot(C, Sc, Tree)

%saveas(figure(1),'/home/mgell/Work/FC/plots/mean_all_FC_clustered_communities_consensus.png');


% % Check for significance by permutation distance between communities of emprical and random nets
% for perm_i = 1:n_perm
%     
%     % Create a degree preserving random graph
%     W0 = null_model_und_sign(W,5,1);
%     
%     % detect communities
%     for part_n = 1:100
%         n0  = size(W0,1);           % number of nodes
%         M0  = 1:n0;                 % initial community affiliations
%         Q0 = -1; Q1 = 0;            % initialize modularity values
%         while Q1-Q0>1e-5;           % while modularity increases
%             Q0 = Q1;                % perform community detection
%             [M0, Q1] = community_louvain(W0, [], M0, 'negative_asym');
%         end
%     end
% 
%     % Calculate partition
%     partition_distance(M, M0)
%     
%     % Save partition dist
%     n_overlap = sum(overlap(overlap == 1))/2;
%     overlaps(perm_i) = n_overlap;
% 
%     clear overlap
% end
% 
% % plot histogram of all permutation overlaps
% figure(1);
% histogram(overlaps,length(unique(overlaps)));
% xticks(unique(overlaps));
% 
% % Calculate empirical overlap wiht prediction edges
% overlap = pred_edges + fing_edges;
% overlap(overlap == 1) = 0;
% overlap(overlap == 2) = 1;
% n_overlap_emp = sum(overlap(overlap == 1))/2;
% 
% % calculate p-value
% p = length(overlaps(overlaps >= n_overlap_emp))/n_perm;



% participation coefficient
% The participation coefficient reflects the extent to which a node is connected to nodes in other modules. 
% https://www.sciencedirect.com/science/article/pii/S1053811916306449#f0010
% In the Guimera and Amaral model, connector hubs are those with high values of participation coefficient and high values of within-module degree z-score.
% Following Guimera and Amaral scheme, we identified connector hubs for each of the three methods as those with simultaneously 
% large values of participation coefficient and within module degree (larger than 0.62 and 1.5, respectively)
[ppos, pneg] = participation_coef_sign(W,M);
Z = module_degree_zscore(W,M,0);

res = [tab array2table(Z) array2table(ppos) array2table(pneg)]

% sort nodes by network
[X,Y,indsort] = grid_communities(M); % can replace M with Sc for consensus clust

% Plot
figure(2);
%set(gcf,'Position',[10 10 730 650])
set(gcf,'Position',[10 10 480 400])
imagesc(W(indsort,indsort));
hold on;
plot(X,Y,'r','linewidth',1);
colorbar();
caxis([-0.1, 0.7]);
%caxis([-0.1, 0.5]); %partial
set(gca,'YTick', [1:length(W)]);
set(gca,'YTickLabel',table2array(tab(indsort,6)));
set(gca,'XTick', [1:length(W)]);
set(gca,'XTickLabel',table2array(tab(indsort,6)));
%ax = gca;
%ax.XAxis.FontSize = 12;
%ax.YAxis.FontSize = 18;
xtickangle(45)
hold off;

%saveas(figure(2),'/home/mgell/Work/FC/plots/mean_all_FC_clustered_communities.png');
%saveas(figure(2),'/home/mgell/Work/FC/plots/mean_all_FC_clustered_communities_21.png');
%writetable(res,'/home/mgell/Work/FC/text_files/res/mean_all_FC_community_prop.csv');
%writetable(res,'/home/mgell/Work/FC/text_files/res/mean_all_FC_community_prop_21.csv');

% Plot consensus
figure(3)
%set(gcf,'Position',[10 10 730 650])
set(gcf,'Position',[10 10 480 400])
imagesc(C(indsort,indsort));
hold on;
colorbar();
caxis([0, 1]);
set(gca,'YTick', [1:length(W)]);
set(gca,'YTickLabel',table2array(tab(indsort,6)));
set(gca,'XTick', [1:length(W)]);
set(gca,'XTickLabel',table2array(tab(indsort,6)));
xtickangle(45)
hold off;

%saveas(figure(3),'/home/mgell/Work/FC/plots/mean_all_FC_clustered_communities_consensus.png');
%saveas(figure(3),'/home/mgell/Work/FC/plots/mean_all_FC_clustered_communities_consensus_21.png');
