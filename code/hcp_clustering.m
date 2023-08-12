
clear

% Perform community detection
n_perm = 1000;

% Packages
addpath(genpath('/home/mgell/Work/Impulsivity_networks/code/'));
addpath(genpath('/home/mgell/Matlab_toolboxes/HierarchicalConsensus-master'));
addpath(genpath('/home/mgell/Matlab_toolboxes/GenLouvain-2.2.0'));

output = '/home/mgell/Work/Impulsivity_networks/hcp/';

% FC & nodes
W = csvread('/home/mgell/Work/Impulsivity_networks/res/replic_5mm_all_mean_s2sFC.csv');

% Nodes
tab = readtable('/home/mgell/Work/Impulsivity_networks/nets/All_localmax.node', 'FileType','text');

% without close neighbours
W([3,17,24],:) = [];
W(:,[3,17,24]) = [];
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
[Sc, Tree] = hierarchicalConsensus(S);
drawHierarchy(Sc, Tree)

C = coclassificationMatrix(S);
consensusPlot(C, Sc, Tree)

%saveas(figure(1),'/home/mgell/Work/FC/hcp_plots/mean_all_FC_clustered_communities_consensus.png');

% Write the order for plotting PET
%writematrix(M,[output 'text_files/res/' nodes 'net_communities.csv'])

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
caxis([-0.2, 0.5]);
%caxis([-0.1, 0.5]); %partial
set(gca,'YTick', [1:length(W)]);
set(gca,'YTickLabel',table2array(tab(indsort,6)));
set(gca,'XTick', [1:length(W)]);
set(gca,'XTickLabel',table2array(tab(indsort,6)));
xtickangle(45)
hold off;

%saveas(figure(2),[output 'hcp_plots/' nodes 'mean_all_FC_clustered_communities.png']);
%saveas(figure(2),[output 'hcp_plots/' nodes 'mean_all_FC_clustered_communities_21.png']);
%saveas(figure(2),[output 'hcp_plots/' nodes 'mean_all_FC_clustered_communities_dorsal_ventral_striatum.png']);
%saveas(figure(2),[output 'hcp_plots/' nodes 'mean_all_FC_clustered_communities_gamma_06.png']);
%writetable(res,[output 'text_files/res/' nodes 'mean_all_FC_community_prop_21.csv']);

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

%saveas(figure(3),[output 'hcp_plots/' 'mean_all_FC_clustered_communities_consensus.png']);
%saveas(figure(3),[output 'hcp_plots/' 'mean_all_FC_clustered_communities_consensus_21.png']);
%saveas(figure(3),[output 'hcp_plots/' 'mean_all_FC_clustered_communities_consensus_dorsal_ventral_striatum.png']);

