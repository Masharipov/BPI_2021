%% ========================================================================
% Masharipov Ruslan, Irina Knyazeva, October, 2021
% Institute of Human Brain of RAS, St. Petersburg, Russia
% Neuroimaging lab
% masharipov@ihb.spb.ru

%% Before running the script:
%  1) Set path to SPM12 (v6906). If you have a later version of SPM,
%     set path to spm_reml.m script from the 'spm_v6906' folder
%  2) Create folder for simulations and copy the 'binary_masks' folder there

%% ========================================================================
close all
clear
% Set path for simulations
sim_path = 'C:\Simulations';

% Load binary masks
load([sim_path filesep 'binary_masks' filesep 'coord.mat']);
act = spm_data_read([sim_path filesep 'binary_masks' filesep 'act_bin.nii'],'xyz',XYZ);
deact = spm_data_read([sim_path filesep 'binary_masks' filesep 'deact_bin.nii'],'xyz',XYZ);
background = spm_data_read([sim_path filesep 'binary_masks' filesep 'background_bin.nii'],'xyz',XYZ);
mask_file = [sim_path filesep 'binary_masks' filesep 'brainmask_bin.nii'];

% Inputs for generation of con-files
iXYZ = cumprod([1,DIM(1:2)'])*XYZ - sum(cumprod(DIM(1:2)')); %coordinates
N =      1000;             % full sample size
PSC =    [0.1 0.2 0.3];    % mean percent signal change (PSC) for (de)activations
SD =     [0.2 0.3 0.4];    % standard deviation of noise
sk =     [0 0.7 -0.7];     % skewness
ku =     [3 2.2 7];        % kurtosis
PSC_bg = 0.045;            % practically non-significant/trivial background  

% Inputs for bayesian estimation
Nstep = [10:10:100 150 200 250 300 350 400 450 500];   % steps for subsampling                       
rep   = 10;                                            % subsampling repetition for each step                          
binary_mask = [sim_path filesep 'binary_masks' filesep 'slice_mask_bin.nii']; % inclusive mask

%% ========================
%    Con-files generation
%  ========================

% Normal distribution (sk = 0, ku = 3)
k = 1; m = 1;
for i = 1:length(PSC)
    for j = 1:length(SD)
        generate_con_files(sim_path,iXYZ,act,deact,background,mask_file,PSC(i),SD(j),sk(k),ku(m),PSC_bg,N,XYZ)
    end
end                 

%% Non-normal distributions

%% (1) weak effect size, high noise (PSC = 0.1, SD = 0.4)
i = 1; j = 3; 
for k = 1:length(sk)
    for m = 1:length(ku)
        if k == 1 & m == 1
            continue
        end
        generate_con_files(sim_path,iXYZ,act,deact,background,mask_file,PSC(i),SD(j),sk(k),ku(m),PSC_bg,N,XYZ)       
    end
end

%% (2) moderate effect size, medium noise (PSC = 0.2, SD = 0.3)
i = 2; j = 2; 
for k = 1:length(sk)
    for m = 1:length(ku)
        if k == 1 & m == 1
            continue
        end
        generate_con_files(sim_path,iXYZ,act,deact,background,mask_file,PSC(i),SD(j),sk(k),ku(m),PSC_bg,N,XYZ)       
    end
end

%% (3) strong effect size, low noise (PSC = 0.3, SD = 0.2)
i = 3; j = 1; 
for k = 1:length(sk)
    for m = 1:length(ku)
        if k == 1 & m == 1
            continue
        end
        generate_con_files(sim_path,iXYZ,act,deact,background,mask_file,PSC(i),SD(j),sk(k),ku(m),PSC_bg,N,XYZ)       
    end
end

%% =========================
%      Con files folders
%  =========================

con_files_all_folders = dir([sim_path filesep 'sim_con_files']);
con_files_all_folders(1:2,:) = [];


%% =======================
%        Estimation
%  =======================

for i = 1:length(con_files_all_folders)
    bayesian_estimation(sim_path,N,Nstep,rep,con_files_all_folders(i).name,binary_mask)
end

%% =======================
%       Ground truth
%  =======================

% Axial slice (z = 36)
slice_mask = spm_data_read(binary_mask,'xyz',XYZ);
% Practically significant effect (Activations, deactivations)
true_act    = act.*slice_mask;        
true_deact  = deact.*slice_mask;      
% Practically non-significant/trivial background (Practical equivalence to the null value)
true_null   = background.*slice_mask; 
% True voxel count for practivally signigicant activations
true_act_count = 100*sum(true_act)/sum(slice_mask);

%% =======================
%         Inference
%  =======================

%% ES threshold = 1 prior SD
for i = 1:length(con_files_all_folders)
    [ROPE_only HDI_ROPE NHST prior_SD] = bayesian_inference(sim_path,Nstep,rep,con_files_all_folders(i).name,binary_mask);
    save([sim_path filesep 'group_stat' filesep con_files_all_folders(i).name '.mat'],'ROPE_only','HDI_ROPE','NHST','prior_SD');
    clear ROPE_only HDI_ROPE NHST prior_SD
end        

%% ES dependencies: 

% normal distribution (sk = 0, ku = 3)
k = 1; m = 1;
% sample size = 200
Nstep = 100; rep = 1;
% ES range
ES = [0:0.001:0.4];

for i = 1:length(PSC)
    for j = 1:length(SD)
        close all
        bayesian_inference_ES_dep(sim_path,Nstep,rep,PSC(i),SD(j),sk(k),ku(m),ES,binary_mask);
    end
end

%% ========================
%  Sample size dependencies
%  ========================

%% Bayesian inference and NHST: Voxel count for 'activated' voxels
LW = 1.5; %line width
scrsz = get(0,'ScreenSize');
close all
for i = 1:length(con_files_all_folders)
    load([sim_path filesep 'group_stat' filesep con_files_all_folders(i).name '.mat']);
    figure(i) 
    set(gcf,'Position',[(scrsz(3)/2-900/2) (scrsz(4)/2-600/2) 400 400])
    errorbar(Nstep,mean(ROPE_only.pos,2),std(ROPE_only.pos,0,2),'-k','MarkerSize',4,'LineWidth',LW); hold on;
    errorbar(Nstep,mean(HDI_ROPE.pos,2),std(HDI_ROPE.pos,0,2),'--k','MarkerSize',4,'LineWidth',LW); hold on;
    errorbar(Nstep,mean(NHST.pos,2),std(NHST.pos,0,2),'-.r','MarkerSize',4,'LineWidth',LW); hold on;
    plot([0 Nstep(end)],[true_act_count true_act_count],'--k','LineWidth',LW);
    xlabel('Sample size'); ylabel('Voxel count'); axis square; axis([0 Nstep(end) 0 mean(NHST.pos(end,:))]);
    title(strrep(con_files_all_folders(i).name, '_', '\_'))
    %legend('ROPE-only','NHST','Location','northwest');
    set(gca,'linewidth',LW)
    print([sim_path filesep 'group_stat' filesep 'Activations_' con_files_all_folders(i).name '.tiff'],'-dtiff','-r600');
    clear ROPE_only HDI_ROPE NHST prior_SD
end

%% Bayesian inference: Voxel count for all voxels
LW = 1.5; %line width
scrsz = get(0,'ScreenSize');
close all
for i = 1:length(con_files_all_folders)
    load([sim_path filesep 'group_stat' filesep con_files_all_folders(i).name '.mat']);
    figure(i)
    set(gcf,'Position',[(scrsz(3)/2-900/2) (scrsz(4)/2-600/2) 400 400])
    errorbar(Nstep,mean(ROPE_only.pos,2),std(ROPE_only.pos,0,2),'-r','MarkerSize',4,'LineWidth',LW); hold on;
    errorbar(Nstep,mean(HDI_ROPE.pos,2),std(HDI_ROPE.pos,0,2),'--r','MarkerSize',4,'LineWidth',LW); hold on;
    errorbar(Nstep,mean(ROPE_only.neg,2),std(ROPE_only.neg,0,2),'-b','MarkerSize',4,'LineWidth',LW); hold on;
    errorbar(Nstep,mean(HDI_ROPE.neg,2),std(HDI_ROPE.neg,0,2),'--b','MarkerSize',4,'LineWidth',LW); hold on;
    errorbar(Nstep,mean(ROPE_only.null,2),std(ROPE_only.null,0,2),'-g','MarkerSize',4,'LineWidth',LW); hold on; 
    errorbar(Nstep,mean(HDI_ROPE.null,2),std(HDI_ROPE.null,0,2),'--g','MarkerSize',4,'LineWidth',LW); hold on; 
    errorbar(Nstep,mean(ROPE_only.lowconf,2),std(ROPE_only.lowconf,0,2),'-k','MarkerSize',4,'LineWidth',LW); hold on; 
    errorbar(Nstep,mean(HDI_ROPE.lowconf,2),std(HDI_ROPE.lowconf,0,2),'--k','MarkerSize',4,'LineWidth',LW); hold on;                
    xlabel('Sample size'); ylabel('Voxel count'); axis square; axis([0 Nstep(end) 0 100]);
    title(strrep(con_files_all_folders(i).name, '_', '\_'))
    %legend('ROPE-only:Act','HDI-ROPE:Act','ROPE-only:Deact','HDI-ROPE:Deact','ROPE-only:Null','HDI-ROPE:Null','ROPE-only:Low-conf','HDI-ROPE:Low-conf'); 
    set(gca,'linewidth',LW)
    print([sim_path filesep 'group_stat' filesep 'Voxel_Count_' con_files_all_folders(i).name '.tiff'],'-dtiff','-r600');
    clear ROPE_only HDI_ROPE NHST prior_SD
end

%% Bayesian inference: Correct decisions (Hit rate or TPR) and Low confidence decisions

LW = 1.5; %line width
scrsz = get(0,'ScreenSize');
close all
for i = 1:length(con_files_all_folders)
    load([sim_path filesep 'group_stat' filesep con_files_all_folders(i).name '.mat']);
    figure(i)
    set(gcf,'Position',[(scrsz(3)/2-900/2) (scrsz(4)/2-600/2) 400 400])
    errorbar(Nstep,mean(ROPE_only.pos_corr,2),std(ROPE_only.pos_corr,0,2),'-r','MarkerSize',4,'LineWidth',LW); hold on;
    errorbar(Nstep,mean(HDI_ROPE.pos_corr,2),std(HDI_ROPE.pos_corr,0,2),'--r','MarkerSize',4,'LineWidth',LW); hold on;
%     errorbar(Nstep,mean(ROPE_only.neg_corr,2),std(ROPE_only.neg_corr,0,2),'-b','MarkerSize',4,'LineWidth',LW); hold on;
%     errorbar(Nstep,mean(HDI_ROPE.neg_corr,2),std(HDI_ROPE.neg_corr,0,2),'--b','MarkerSize',4,'LineWidth',LW); hold on;
    errorbar(Nstep,mean(ROPE_only.null_corr,2),std(ROPE_only.null_corr,0,2),'-g','MarkerSize',4,'LineWidth',LW); hold on; 
    errorbar(Nstep,mean(HDI_ROPE.null_corr,2),std(HDI_ROPE.null_corr,0,2),'--g','MarkerSize',4,'LineWidth',LW); hold on; 
    errorbar(Nstep,mean(ROPE_only.lowconf,2),std(ROPE_only.lowconf,0,2),'-k','MarkerSize',4,'LineWidth',LW); hold on; 
    errorbar(Nstep,mean(HDI_ROPE.lowconf,2),std(HDI_ROPE.lowconf,0,2),'--k','MarkerSize',4,'LineWidth',LW); hold on;                
    xlabel('Sample size'); ylabel('Correct and Low confidence decisions'); axis square; axis([0 Nstep(end) 0 100]);
    title(strrep(con_files_all_folders(i).name, '_', '\_'))
    %legend('ROPE-only:Act','HDI-ROPE:Act','ROPE-only:Deact','HDI-ROPE:Deact','ROPE-only:Null','HDI-ROPE:Null','ROPE-only:Low-conf','HDI-ROPE:Low-conf'); 
    set(gca,'linewidth',LW)
    print([sim_path filesep 'group_stat' filesep 'Corr_dec_' con_files_all_folders(i).name '.tiff'],'-dtiff','-r600');
    clear ROPE_only HDI_ROPE NHST prior_SD    
end

%% Bayesian inference: Incorrect decisions (False alarm or FPR)

LW = 1.5; %line width
scrsz = get(0,'ScreenSize');
close all
for i = 1:length(con_files_all_folders)
    load([sim_path filesep 'group_stat' filesep con_files_all_folders(i).name '.mat']);
    figure(i)
    errorbar(Nstep,mean(ROPE_only.pos_incorr,2),std(ROPE_only.pos_incorr,0,2),'-r','MarkerSize',4,'LineWidth',LW); hold on;
    errorbar(Nstep,mean(HDI_ROPE.pos_incorr,2),std(HDI_ROPE.pos_incorr,0,2),'--r','MarkerSize',4,'LineWidth',LW); hold on;
    errorbar(Nstep,mean(ROPE_only.neg_incorr,2),std(ROPE_only.neg_incorr,0,2),'-b','MarkerSize',4,'LineWidth',LW); hold on;
    errorbar(Nstep,mean(HDI_ROPE.neg_incorr,2),std(HDI_ROPE.neg_incorr,0,2),'--b','MarkerSize',4,'LineWidth',LW); hold on;
    errorbar(Nstep,mean(ROPE_only.null_incorr,2),std(ROPE_only.null_incorr,0,2),'-g','MarkerSize',4,'LineWidth',LW); hold on; 
    errorbar(Nstep,mean(HDI_ROPE.null_incorr,2),std(HDI_ROPE.null_incorr,0,2),'--g','MarkerSize',4,'LineWidth',LW); hold on;            
    xlabel('Sample size'); ylabel('Incorrect decisions'); axis square; axis([0 Nstep(end) 0 5]);
    title(strrep(con_files_all_folders(i).name, '_', '\_'))
    %legend('ROPE-only:Act','HDI-ROPE:Act','ROPE-only:Null','ROPE-only:Deact','HDI-ROPE:Deact','HDI-ROPE:Null','ROPE-only:Low-conf','HDI-ROPE:Low-conf');
    set(gca,'linewidth',LW)
    print([sim_path filesep 'group_stat' filesep 'Incorr_dec_' con_files_all_folders(i).name '.tiff'],'-dtiff','-r600');
    clear ROPE_only HDI_ROPE NHST prior_SD    
end
   
