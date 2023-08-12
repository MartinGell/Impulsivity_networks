% gradient calc

clear

addpath(genpath('/home/mgell/Matlab_toolboxes/BrainSpace/'))
%addpath(genpath('/data/project/impulsivity/toolboxes/BrainSpace/'))

tab = readtable('/home/mgell/Work/Impulsivity_networks/nets/All_localmax.node', 'FileType','text');
tab([2,17,24],:) = [];

% similarity
s2v_similarity = csvread('/home/mgell/Work/Impulsivity_networks/res/s2v_corr.csv');
%s2v_similarity = csvread('/home/mgell/Work/FC/text_files/res/s2v_similarity.csv');
%s2v_similarity = csvread('/data/project/impulsivity/FC/s2v_similarity.csv');
%s2v_similarity = csvread('/data/project/impulsivity/FC/s2v_corr.csv');

s2v_similarity(16,:) = []; s2v_similarity(:,16) = [];

%s2v_similarity = s2v_similarity - diag(diag(s2v_similarity));

figure(1);
set(gcf,'Position',[10 10 730 650])
imagesc(s2v_similarity);
hold on;
colorbar();
%caxis([0.8, 1]);
set(gca,'YTick', [1:22]);
set(gca,'YTickLabel',table2array(tab(:,6)));
set(gca,'XTick', [1:22]);
set(gca,'XTickLabel',table2array(tab(:,6)));
xtickangle(45)
hold off;

% Construct the gradients
gm = GradientMaps('kernel','na','approach','pca','n_components',10);
gm = gm.fit(s2v_similarity, 'sparsity', 70);

gr = gm.gradients{1};
%writematrix(gr,'/home/mgell/Work/FC/text_files/res/grawritematrix(gr,'/home/mgell/Work/FC/text_files/res/gradient_values_pca_70_zscore.csv');dient_values_pca_80_zscore.csv');
%writematrix(gr,'/home/mgell/Work/FC/text_files/res/gradient_values_pca_70_zscore.csv');

scree_plot(gm.lambda{1});
%saveas(figure(2),'/data/project/impulsivity/FC/gradients_res.png');
%saveas(figure(2),'/home/mgell/Work/FC/plots/pca_gradients_lambda_res_70.png');

% 2D
%gradient_in_euclidean(gm.gradients{1}(:,1:2));
%saveas(figure(2),'/data/project/impulsivity/FC/g1_by_g2.png');
%saveas(figure(1),'/home/mgell/Work/test/g1_by_g2_dm.png');

% 3D
gradient_in_euclidean(gm.gradients{1}(:,1:3));
%saveas(figure(2),'/home/mgell/Work/FC/plots/dm_all_gradients3D.png');
% 
% xxx = rand(22,1)*0.01;
% yyy = rand(22,1)*0.02;
% 
% figure(3); scatter(gm.gradients{1}(:,2),gm.gradients{1}(:,3)); 
% text(gm.gradients{1}(:,2)+xxx,gm.gradients{1}(:,3)+xxx,table2array(tab(:,5)));
% text(gm.gradients{1}(:,1),gm.gradients{1}(:,2),table2array(tab(:,5)));
% 
% gr = gm.gradients{1};
% save('/data/project/impulsivity/FC/gradients_res.mat', 'gr');
% writematrix(gr,'/home/mgell/Work/FC/text_files/res/gradient_values.csv');
% 
% lb = gm.lambda{1};
% save('/data/project/impulsivity/FC/gradients_res_lambda.mat', 'lb');
% 
% gm
