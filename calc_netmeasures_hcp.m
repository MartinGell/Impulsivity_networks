
clear


% Assemble seed to voxel FC mats from all subjects and average -> create a
% nifti with average seed2voxel FC.
% Using mat files becuase readtable and readmatrix are very slow reading in
% the .csv files. Possibly concattenating first would speed this up.
% readcsv just loads it in as one long vector.


sample = 'all'       % always 'all' for HCP
nodes = 'replic_5mm_'% replic_5mm_ or replic_6mm_ or extended_5mm_
mat = 1              % input as matrix or vector? With current pipeline should be 1
denoise = 0          % regress confounds or not?
partial = 0          % average normal or partial cor?



% -------------------- SET UP -----------------------
% paths
cd('/home/mgell/Work/To_juseless/hcp_analysis');
addpath(genpath('/home/mgell/Work/To_juseless/hcp_analysis'));
addpath(genpath('/home/mgell/Work/FC/code'));

nii_toolbox = '/home/mgell/Matlab_toolboxes/NIFTI_toolbox-master/';

M = readmatrix('/home/mgell/Work/FC/hcp/text_files/res/replic_5mm_net_communities.csv');

% subject FC location and table with nodes and names
switch nodes
    case 'replic_5mm_'
        inputdir = '/home/mgell/Work/To_juseless/Nevena/hcp_fingerprint/output/hcp/All_localmax/';
        tab = readtable('/home/mgell/Work/FC/nets/All_localmax.node', 'FileType','text');
        tab([2,24],:) = [];
    case 'replic_6mm_'
        inputdir = '/home/mgell/Work/To_juseless/hcp_fingerprint/output/hcp/imp_all_net_RL/';
        tab = readtable('/home/mgell/Work/FC/nets/All_localmax.node', 'FileType','text');
    case 'extended_5mm_'
        inputdir = '/home/mgell/Work/To_juseless/hcp_fingerprint/output/hcp/Extended_net_localmax/';
        tab = readtable('/home/mgell/Work/FC/nets/All_localmax_brainnet_extended.node', 'FileType','text');
end

% output directory
outputdir = fullfile('/home/mgell/Work/FC/hcp/text_files/res/'); 

if denoise == 1
    % confounds to regress
    BV = readtable('/home/mgell/Work/FC/text_files/brainsize.csv');
    BV = table2array(BV(1,:));
    subs = readtable('/home/mgell/Work/FC/text_files/final_subs.csv');
    subs.BV = BV';
    subs(subs.age < 19,:) = [];
    subs = subs(1:100,:);
    age = table2array(subs(:,5));
    FD = table2array(subs(:,9));
    conf = [FD table2array(subs(:,17)) age];
end
% ----------------- END OF SET UP --------------------



% Load sample to average over
%data = dir([inputdir '*.mat']);
sample_subs = table2array(readtable(fullfile(pwd,'samples', [sample '.csv'])));
%sample_subs = sample_subs(1:100,1);

if exist('all_data','var')
    clear('all_data')    
end

i = 1;

for file_i = 1:length(sample_subs)
    % load subject
    sub_i = sample_subs(file_i);
    sub_i = [num2str(sub_i) '_conn_s5mm_filt0_010_1_confWM_CSF_GS_confdrv_maskCAT12_02_dt11'];
    try
    subFC = load(fullfile(inputdir, [sub_i '.mat']));
    catch
        fprintf('Subject %s missing \n', sub_i)
        continue
    end
    
    if partial == 0
        subFC = subFC.mean_mat_c;
    else
        subFC = subFC.meanpartial_mat_c;
    end
    
    % make sure all subjects have the correct lenght in second dimension (n of vox)
%     [seed_N,vox_N] = size(subFC);
%     if vox_N ~= 156605
%         fprintf('subject %s has different number of voxels than group', file.name)
%         continue
%     elseif seed_N ~= 17
%         fprintf('subject %s has different number of seeds than group', file.name)
%         continue
%     end

    % check for NaNs
    N = sum(isnan(subFC(:)));
    if N ~= 0
        fprintf('subject %s has %i NaNs removing... \n', sub_i{1}, N)
        continue
    end

    % save
    W = atanh(subFC);
    W([2,24],:) = [];
    W(:,[2,24]) = [];
    W(find(eye(length(W)))) = zeros(1,length(22)); % replace diag with 1
    
    [X,Y,indsort] = grid_communities(M);
    %W = W(indsort,indsort);
        
    mean_net = ones(4,4);
    for net_i = 1:4
        nodes_i = indsort(M(:,1) == net_i); % WHEN net_i == net_j have to only take upper triangle otherwise fc will be approx zero 
        
        for net_j = 1:4
            nodes_j = indsort(M(:,1) == net_j);
            B = W([nodes_i],[nodes_j]);
            if net_i == net_j
                B = B(triu(true(size(B)),+1))';
            end
            mean_net(net_i,net_j) = mean(B,"all");
            
            clear B
        end
    end
  
    all_data(:,i) = mean_net(triu(true(size(mean_net))))';

% 
%     [ppos, pneg] = participation_coef_sign(W,M);    
%     [C_pos,C_neg,Ctot_pos,Ctot_neg] = clustering_coef_wu_sign(W);
%     
%     clust_all(i) = Ctot_pos; 
%     eff_all(i) = efficiency_wei(W);
%     Z_all(:,i) = module_degree_zscore(W,M,0);
%     ppos_all(:,i) = ppos;

    % update counter
    i = i+1;
end

% save resulting mean seed2vox FC
%dlmwrite(strjoin({saveit atlas '_mean_fc.txt'},''), mean_fc, 'delimiter', ' '); 

all_data = vertcat(sample_subs',all_data);

if partial == 0
    %writematrix(Z_all', fullfile(outputdir, [nodes sample '_Z.csv']));
    %writematrix(ppos_all', fullfile(outputdir, [nodes sample '_ppos.csv']));
    %writematrix(clust_all', fullfile(outputdir, [nodes sample '_clust.csv']));
    %writematrix(eff_all', fullfile(outputdir, [nodes sample '_eff.csv']));
    writematrix(all_data', fullfile(outputdir, [nodes sample '_mean_netFC.csv']));
else
    writematrix(mean_s2sFC, fullfile(outputdir, [nodes sample '_partial_mean_s2sFC.csv']));
end
%mean_s2sFC = mean_s2sFC - diag(diag(mean_s2sFC));




% -------------------- PLOT -----------------------
figure(1);
set(gcf,'Position',[10 10 730 650])
imagesc(W);
hold on;
colorbar();
caxis([-0.2, 0.5]);
%caxis([-0.1, 0.4]); %partial
set(gca,'YTick', [1:length(W)]);
set(gca,'YTickLabel',table2array(tab(:,6)));
set(gca,'XTick', [1:length(W)]);
set(gca,'XTickLabel',table2array(tab(:,6)));
xtickangle(45)
hold off;

if partial == 0
    saveas(figure(1), fullfile(outputdir, [nodes sample '_plot_mean_s2sFC.png']));
else
    saveas(figure(1), fullfile(outputdir, [nodes sample '_plot_partial_mean_s2sFC.png']));
end

