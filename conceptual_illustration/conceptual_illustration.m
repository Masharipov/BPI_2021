%% ========================================================================
% Masharipov Ruslan, October, 2021
% Institute of Human Brain of RAS, St. Petersburg, Russia
% Neuroimaging lab
% masharipov@ihb.spb.ru

%% Before running the script:
%  Set path to SPM12 (v6906). If you have a later version of SPM,
%  set path to spm_reml.m script from the 'spm_v6906' folder

%% ========================================================================
close all
clear

% set path for simulations
sim_path = 'C:\conceptual_illustration';

status = exist([sim_path filesep 'illustation']);
if status ~= 7
    mkdir([sim_path filesep 'illustation']); 
else
    rmdir([sim_path filesep 'illustation'],'s');
    mkdir([sim_path filesep 'illustation']); 
end

% activation areas (number '1')
activ = zeros(50,50);   
activ(8:42,15:17) = 1;
activ(13:15,12:17) = 1;
activ(39:42,12:21) = 1;

% null areas (number '0')
null = zeros(50,50);
null(11:39,31:34) = 1;
null(11:39,39:42) = 1;
null(8:11,32:41) = 1;
null(39:42,32:41) = 1;

% low confidence areas
low_conf = ones(50,50) - activ - null;   

% ground truth
ground_truth_mean = activ*0.1 + null*0 + low_conf*0.01;
ground_truth_sd   = activ*0.37 + null*0.06 + low_conf*0.37;
figure(1)
subplot(1,2,1); imagesc(ground_truth_mean); title('Ground truth mean'); axis square; colorbar;
subplot(1,2,2); imagesc(ground_truth_sd);   title('Ground truth SD');   axis square; colorbar;

% add noise, smooth and generate *.nii files
N = 100; % sample size
for i=1:N 
    %noise
    noise = normrnd(zeros(50,50),ground_truth_sd,50,50); 
    %image
    image = ground_truth_mean + noise;
    s_image(:,:,i) = imGaussFilter(image,2);
    nii_image = make_nii(s_image(:,:,i));
    save_nii(nii_image,[sim_path  filesep 'illustation' filesep  num2str(i,'%03.f') '_image.nii']); 
    clear noise image nii_image
end

mean_image   = mean(s_image,3);
sd_image     = std(s_image,0,3);
se_image     = sd_image/sqrt(N);
t_image      = mean_image./se_image;

%% Bayesian estimation
mkdir([sim_path filesep 'illustation' filesep 'one_sample'])
matlabbatch{1}.spm.stats.factorial_design.dir = {[sim_path  filesep 'illustation' filesep 'one_sample']};
for j=1:N
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{j,1} = [sim_path filesep 'illustation' filesep num2str(j,'%03.f') '_image.nii,1'];   
end
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{3}.spm.stats.fmri_est.method.Bayesian2 = 1;
spm_jobman('run',matlabbatch);

%% Bayesian inference
cd([sim_path filesep 'illustation' filesep 'one_sample'])
load('SPM.mat')

% read posterior beta
XYZ  = SPM.xVol.XYZ;
iXYZ = cumprod([1,SPM.xVol.DIM(1:2)'])*XYZ - sum(cumprod(SPM.xVol.DIM(1:2)'));
cB    = spm_data_read(SPM(1).VCbeta,'xyz',XYZ); 
c = 1;
cB = c*cB;

% compute posterior variance
VcB   = c'*SPM.PPM.Cby*c;
for j = 1:length(SPM.PPM.l)
    l   = spm_data_read(SPM.VHp(j),'xyz',XYZ);
    VcB = VcB + (c'*SPM.PPM.dC{j}*c)*(l - SPM.PPM.l(j));
end             

% prior SD
prior_SD = full(sqrt(c'*SPM.PPM.Cb*c));

% effect size threshold
ES = prior_SD;

% mask
mask = spm_read_vols(spm_vol('mask.nii'));

% log posterior odds (LPOs)
PostOdds_null = (normcdf(ES,cB,sqrt(VcB)) - normcdf(-ES,cB,sqrt(VcB)))...
            ./(1 - normcdf(ES,cB,sqrt(VcB)) + normcdf(-ES,cB,sqrt(VcB)));
              
LPO_null = log(PostOdds_null); % evidence for H-null
LPO_alt = -LPO_null;     	   % evidence for H-alt

LPO_alt_img = mask;
LPO_alt_img(iXYZ) = LPO_alt;

figure(2)
subplot(1,2,1); imagesc(t_image);     title('Classical NSHT'); axis square; 
caxis([-2.5 6]); h1 = colorbar; h1.Title.String = 'T-value';
subplot(1,2,2); imagesc(LPO_alt_img); title('Bayesian inference'); axis square; 
caxis([-3.5 3.5]); h2 = colorbar; h2.Title.String = 'LPO';

cd(sim_path)



