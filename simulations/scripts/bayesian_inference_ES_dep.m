function bayesian_inference_ES_dep(sim_path,Nstep,rep,PSC,SD,sk,ku,ES,binary_mask)

tic
load([sim_path filesep 'group_stat' filesep 'PSC_(' num2str(PSC) ')_SD_(' num2str(SD) ')_Sk_(' num2str(sk) ')_Ku_(' num2str(ku) ')' ...
     filesep 'Sample_' num2str(Nstep,'%03d') '_Rep_' num2str(rep,'%03d') filesep 'SPM.mat']);
cd(SPM.swd)
XYZ = SPM.xVol.XYZ;    % coordites
c = 1;                 % contrast    
% Load betas
class_B = c*spm_data_read(SPM.Vbeta,'xyz',XYZ);      % classical 
cB      = c*spm_data_read(SPM.VCbeta,'xyz',XYZ);     % bayesian        
% Load variance
class_l   = spm_data_read(SPM.VResMS,'xyz',XYZ);     % get hyperparamters
class_Vc  = c'*SPM.xX.Bcov*c;
SE  = sqrt(class_l*class_Vc);                        % and standard error
Z   = class_B./SE;                
% Compute posterior variance
VcB   = c'*SPM.PPM.Cby*c;
for k = 1:length(SPM.PPM.l)
    l   = spm_data_read(SPM.VHp(k),'xyz',XYZ);              % hyperparameter
    VcB = VcB + (c'*SPM.PPM.dC{k}*c)*(l - SPM.PPM.l(k));    % Taylor approximation
end  

prior_SD = full(sqrt(SPM.PPM.Cb));
        
% 95% HDI 
HDImax = spm_invNcdf(0.975,cB,VcB);
HDImin = spm_invNcdf(0.025,cB,VcB); 

% Ground truth
total = length(cB); 
slice_mask  = spm_data_read(binary_mask,'xyz',XYZ);
act         = spm_data_read([sim_path filesep 'binary_masks' filesep 'act_bin.nii'],'xyz',XYZ);
deact       = spm_data_read([sim_path filesep 'binary_masks' filesep 'deact_bin.nii'],'xyz',XYZ);
background  = spm_data_read([sim_path filesep 'binary_masks' filesep 'background_bin.nii'],'xyz',XYZ);
% Practically significant effect (Activations, deactivations)
true_act    = act.*slice_mask;        
true_deact  = deact.*slice_mask;      
% Practically non-significant/trivial background (Practical equivalence to the null value)
true_null   = background.*slice_mask; 

% ES threshold (ROPE radius) dependency
for thr = 1:length(ES)
    
    ROPE_max = ES(thr);
    ROPE_min = -ES(thr);
        
    % Posterior probability
    post_pos = (normcdf(-ES(thr),-cB,sqrt(VcB)));
    post_neg = (normcdf(-ES(thr),cB,sqrt(VcB)));
    post_null = (normcdf(ES(thr),cB,sqrt(VcB)) - normcdf(-ES(thr),cB,sqrt(VcB)));
    
    % Voxel count            
    ROPE_only(1).pos(thr)  = 100*sum(post_pos>=0.95)/total;
    ROPE_only(1).neg(thr)  = 100*sum(post_neg>=0.95)/total;
    ROPE_only(1).null(thr) = 100*sum(post_null>=0.95)/total;
    ROPE_only(1).lowconf(thr) = 100 - ROPE_only(1).pos(thr) - ROPE_only(1).neg(thr) - ROPE_only(1).null(thr);
    HDI_ROPE(1).pos(thr)  = 100*sum(HDImin>=ROPE_max)/total;
    HDI_ROPE(1).neg(thr)  = 100*sum(HDImax<=ROPE_min)/total;
    HDI_ROPE(1).null(thr) = 100*sum(HDImin>=ROPE_min & HDImax<=ROPE_max)/total;
    HDI_ROPE(1).lowconf(thr) = 100 - HDI_ROPE(1).pos(thr) - HDI_ROPE(1).neg(thr) - HDI_ROPE(1).null(thr);

    % Correct decisions (Hit rate or TPR)
    ROPE_only(1).pos_corr(thr)  = 100*sum((post_pos>=0.95).*true_act)/sum(true_act);
    ROPE_only(1).neg_corr(thr)  = 100*sum((post_neg>=0.95).*true_deact)/sum(true_deact);
    ROPE_only(1).null_corr(thr) = 100*sum((post_null>=0.95).*true_null)/sum(true_null); 
    HDI_ROPE(1).pos_corr(thr)   = 100*sum((HDImin>=ROPE_max).*true_act)/sum(true_act);
    HDI_ROPE(1).neg_corr(thr)   = 100*sum((HDImax<=ROPE_min).*true_deact)/sum(true_deact);
    HDI_ROPE(1).null_corr(thr)  = 100*sum((HDImin>=ROPE_min & HDImax<=ROPE_max).*true_null)/sum(true_null);
        
    % Incorrect decisions (False alarm or FPR)
    ROPE_only_incorr_pos  = sum(((post_pos>=0.95)  - true_act)==1);     
    ROPE_only_incorr_neg  = sum(((post_neg>=0.95)  - true_deact)==1);   
    ROPE_only_incorr_null = sum(((post_null>=0.95) - true_null)==1);    
    HDI_ROPE_incorr_pos   = sum(((HDImin>=ROPE_max)  - true_act)==1);   
    HDI_ROPE_incorr_neg   = sum(((HDImax<=ROPE_min)  - true_deact)==1); 
    HDI_ROPE_incorr_null  = sum(((HDImin>=ROPE_min & HDImax<=ROPE_max) - true_null)==1);

    ROPE_only(1).pos_incorr(thr)  = 100*ROPE_only_incorr_pos/(total - sum(true_act));
    ROPE_only(1).neg_incorr(thr)  = 100*ROPE_only_incorr_neg/(total - sum(true_deact));
    ROPE_only(1).null_incorr(thr) = 100*ROPE_only_incorr_null/(total - sum(true_null)); 
    HDI_ROPE(1).pos_incorr(thr)   = 100*HDI_ROPE_incorr_pos/(total - sum(true_act));
    HDI_ROPE(1).neg_incorr(thr)   = 100*HDI_ROPE_incorr_neg/(total - sum(true_deact)); 
    HDI_ROPE(1).null_incorr(thr)  = 100*HDI_ROPE_incorr_null/(total - sum(true_null)); 
end


%% Correct (Hit rate or TPR) and incorrect decisions (False alarm or FPR)
LW = 2; %line width
scrsz = get(0,'ScreenSize');
figure
set(gcf,'Position',[(scrsz(3)/2-900/2) (scrsz(4)/2-600/2) 400 400])
plot(ES,ROPE_only.pos_corr,'-r','MarkerSize',4,'LineWidth',LW); hold on;
plot(ES,ROPE_only.pos_incorr,'-.r','MarkerSize',4,'LineWidth',LW); hold on;
plot(ES,ROPE_only.null_corr,'-g','MarkerSize',4,'LineWidth',LW); hold on; 
plot(ES,ROPE_only.null_incorr,'-.g','MarkerSize',4,'LineWidth',LW); hold on; 
plot([prior_SD prior_SD],[0 100],'--k','LineWidth',LW);       
xlabel('ES threshold,%'); ylabel('Correct and incorrect decisions'); axis square; axis([0 ES(end) 0 100]);
set(gca,'linewidth',1.5)
print([sim_path filesep 'group_stat' filesep 'PSC_(' num2str(PSC) ')_SD_(' num2str(SD) ')_Sk_(' num2str(sk) ')_Ku_(' num2str(ku) ')' ...
       '_N' num2str(Nstep) '_ES_dependency.tiff'],'-dtiff','-r600');

toc

end
