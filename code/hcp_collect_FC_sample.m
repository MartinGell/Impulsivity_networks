
clear


% Assemble seed to voxel FC mats from all subjects and average -> create a
% nifti with average seed2voxel FC.
% Using mat files becuase readtable and readmatrix are very slow reading in
% the .csv files. Possibly concattenating first would speed this up.
% readcsv just loads it in as one long vector.


sample = 'all'       % always 'all' for HCP
mat = 1              % input as matrix or vector? With current pipeline should be 1
denoise = 0          % regress confounds or not?
partial = 0          % average normal or partial cor?



% -------------------- SET UP -----------------------
% paths
cd('/home/mgell/Work/Impulsivity_networks/code/');
addpath(genpath('/home/mgell/Work/Impulsivity_networks/code/'));

nii_toolbox = '/home/mgell/Matlab_toolboxes/NIFTI_toolbox-master/';

% subject FC location
inputdir = '/home/mgell/Work/Impulsivity_networks/input/';

% output directory
outputdir = '/home/mgell/Work/Impulsivity_networks/res/'; 

% table with nodes and names
tab = readtable('/home/mgell/Work/Impulsivity_networks/nets/All_localmax.node', 'FileType','text');

if denoise == 1
    % confounds to regress
    BV = readtable('/home/mgell/Work/FC/text_files/brainsize.csv');  % needs a csv with varaibles to denoise*N with subject id as column header
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
% list of subjects to average FC over
sample_subs = table2array(readtable(fullfile('/home/mgell/Work/Impulsivity_networks/Impulsivity_networks/text_files', [sample '.csv']))); 


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
    
    % save
    if mat == 1
        all_data(:,i) = atanh(subFC(triu(true(size(subFC)),+1))');
    else
        all_data(:,i) = atanh(subFC);
    end

    % check for NaNs
    N = sum(isnan(subFC(:)));
    if N ~= 0
        fprintf('subject %s has %i NaNs removing... \n', sub_i{1}, N)
        continue
    end
    
    % update counter
    i = i+1;
end

% average all seed-vox corrs
mean_s2sFC = mean(all_data,2);
mean(all_data(1,:)) == mean_s2sFC(1,1) % check mean is correctly applied

% regress confounds from FC
if denoise == 1
    % detrend confounds
    for rc=1:size(conf,2)
       try
           conf(:,rc) = detrend(conf(:,rc),1);
       catch MEME
       end
    end

    all_data = all_data';
    
    conf = [conf, ones(size(conf,1),1)];
    bb = conf\all_data;
    deconf_data = all_data-conf*bb;
    
    % average all seed-vox corrs
    Xdeconf_mean_s2sFC = mean(deconf_data,1);
    mean(deconf_data(:,1)) == Xdeconf_mean_s2sFC(1,1) % check mean is correctly applied
    
    % put back mean
    new_mean_s2sFC = Xdeconf_mean_s2sFC' + bb(4,:)'; % should this be intercept or mean > mean_s2sFC;
    mean_s2sFC = new_mean_s2sFC;
end

% reasemble FC vecs to mat
%mean_s2sFC = tanh(mean_s2sFC);
FC_mat = reasemble_FC_mat(mean_s2sFC,size(tab,1),1); % make usure to use the same tril/triu here and in function
mean_s2sFC = FC_mat;

% save resulting mean seed2vox FC
%dlmwrite(strjoin({saveit atlas '_mean_fc.txt'},''), mean_fc, 'delimiter', ' '); 

if partial == 0
    writematrix(mean_s2sFC, fullfile(outputdir, [sample '_mean_s2sFC.csv']));
else
    writematrix(mean_s2sFC, fullfile(outputdir, [sample '_partial_mean_s2sFC.csv']));
end
%mean_s2sFC = mean_s2sFC - diag(diag(mean_s2sFC));




% -------------------- PLOT -----------------------
figure(1);
set(gcf,'Position',[10 10 730 650])
imagesc(mean_s2sFC);
hold on;
colorbar();
caxis([-0.2, 0.5]);
%caxis([-0.1, 0.4]); %partial
set(gca,'YTick', [1:length(subFC)]);
set(gca,'YTickLabel',table2array(tab(:,6)));
set(gca,'XTick', [1:length(subFC)]);
set(gca,'XTickLabel',table2array(tab(:,6)));
xtickangle(45)
hold off;

if partial == 0
    saveas(figure(1), fullfile(outputdir, [sample '_plot_mean_s2sFC.png']));
else
    saveas(figure(1), fullfile(outputdir, [sample '_plot_partial_mean_s2sFC.png']));
end

