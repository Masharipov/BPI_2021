function [ROPE_only HDI_ROPE NHST prior_SD] = bayesian_inference(sim_path,Nstep,rep,con_files_folder,binary_mask)

for i = 1:length(Nstep) 
    for j = 1:rep
        load([sim_path filesep 'group_stat' filesep con_files_folder filesep 'Sample_' num2str(Nstep(i),'%03d') '_Rep_' num2str(j,'%03d') filesep 'SPM.mat']);
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
        % ES threshold (ROPE radius)
        ES = full(sqrt(SPM.PPM.Cb));        % prior SD
        prior_SD(i,j) = ES;
        ROPE_max = ES;
        ROPE_min = -ES;       
        % 95% HDI 
        HDImax = spm_invNcdf(0.975,cB,VcB);
        HDImin = spm_invNcdf(0.025,cB,VcB);             
        % Posterior probability
        post_pos = (normcdf(-ES,-cB,sqrt(VcB)));
        post_neg = (normcdf(-ES,cB,sqrt(VcB)));
        post_null = (normcdf(ES,cB,sqrt(VcB)) - normcdf(-ES,cB,sqrt(VcB)));
        
        % ROPE-only (PP>0.95 or LPO>3)        
        ROPE_only(1).statmap(i).sample(j).rep.pos  = (post_pos>=0.95);
        ROPE_only(1).statmap(i).sample(j).rep.neg  = (post_neg>=0.95);
        ROPE_only(1).statmap(i).sample(j).rep.null = (post_null>=0.95);
                       
        % HDI-ROPE (overlap between 95% HDI and ROPE)      
        HDI_ROPE(1).statmap(i).sample(j).rep.pos  = (HDImin>=ROPE_max);
        HDI_ROPE(1).statmap(i).sample(j).rep.neg  = (HDImax<=ROPE_min);
        HDI_ROPE(1).statmap(i).sample(j).rep.null = (HDImin>=ROPE_min & HDImax<=ROPE_max);
          
        % Classical NHST voxel-wise pFWE<0.05
        df = [1 SPM.xX.erdf];
        S    = SPM.xVol.S;                  %-search Volume {voxels}
        R    = SPM.xVol.R;                  %-search Volume {resels}
        u = spm_uc(0.025,df,'T',R,1,S);
        NHST(1).statmap(i).sample(j).rep.pos = (Z>=u);
        NHST(1).statmap(i).sample(j).rep.neg = (Z<=-u);
        
        % Voxel count
        total = length(cB);              
        ROPE_only(1).pos(i,j)  = 100*sum(post_pos>=0.95)/total;
        ROPE_only(1).neg(i,j)  = 100*sum(post_neg>=0.95)/total;
        ROPE_only(1).null(i,j) = 100*sum(post_null>=0.95)/total;
        ROPE_only(1).lowconf(i,j) = 100 - ROPE_only(1).pos(i,j) - ROPE_only(1).neg(i,j) - ROPE_only(1).null(i,j);
        HDI_ROPE(1).pos(i,j)  = 100*sum(HDImin>=ROPE_max)/total;
        HDI_ROPE(1).neg(i,j)  = 100*sum(HDImax<=ROPE_min)/total;
        HDI_ROPE(1).null(i,j) = 100*sum(HDImin>=ROPE_min & HDImax<=ROPE_max)/total;
        HDI_ROPE(1).lowconf(i,j) = 100 - HDI_ROPE(1).pos(i,j) - HDI_ROPE(1).neg(i,j) - HDI_ROPE(1).null(i,j);
        NHST(1).pos(i,j) = 100*sum(Z>=u)/total;
        NHST(1).neg(i,j) = 100*sum(Z<=-u)/total;
        
        % Ground truth
        slice_mask  = spm_data_read(binary_mask,'xyz',XYZ);
        act         = spm_data_read([sim_path filesep 'binary_masks' filesep 'act_bin.nii'],'xyz',XYZ);
        deact       = spm_data_read([sim_path filesep 'binary_masks' filesep 'deact_bin.nii'],'xyz',XYZ);
        background  = spm_data_read([sim_path filesep 'binary_masks' filesep 'background_bin.nii'],'xyz',XYZ);
        % Practically significant effect (Activations, deactivations)
        true_act    = act.*slice_mask;        
        true_deact  = deact.*slice_mask;      
        % Practically non-significant/trivial background (Practical equivalence to the null value)
        true_null   = background.*slice_mask; 
        
        % Correct decisions (Hit rate or TPR)
        ROPE_only(1).pos_corr(i,j)  = 100*sum((post_pos>=0.95).*true_act)/sum(true_act);
        ROPE_only(1).neg_corr(i,j)  = 100*sum((post_neg>=0.95).*true_deact)/sum(true_deact);
        ROPE_only(1).null_corr(i,j) = 100*sum((post_null>=0.95).*true_null)/sum(true_null); 
        HDI_ROPE(1).pos_corr(i,j)   = 100*sum((HDImin>=ROPE_max).*true_act)/sum(true_act);
        HDI_ROPE(1).neg_corr(i,j)   = 100*sum((HDImax<=ROPE_min).*true_deact)/sum(true_deact);
        HDI_ROPE(1).null_corr(i,j)  = 100*sum((HDImin>=ROPE_min & HDImax<=ROPE_max).*true_null)/sum(true_null);
        
        % Incorrect decisions (False alarm or FPR)
        ROPE_only_incorr_pos  = sum(((post_pos>=0.95)  - true_act)==1);     
        ROPE_only_incorr_neg  = sum(((post_neg>=0.95)  - true_deact)==1);   
        ROPE_only_incorr_null = sum(((post_null>=0.95) - true_null)==1);    
        HDI_ROPE_incorr_pos   = sum(((HDImin>=ROPE_max)  - true_act)==1);   
        HDI_ROPE_incorr_neg   = sum(((HDImax<=ROPE_min)  - true_deact)==1); 
        HDI_ROPE_incorr_null  = sum(((HDImin>=ROPE_min & HDImax<=ROPE_max) - true_null)==1);
        
        ROPE_only(1).pos_incorr(i,j)  = 100*ROPE_only_incorr_pos/(total - sum(true_act));
        ROPE_only(1).neg_incorr(i,j)  = 100*ROPE_only_incorr_neg/(total - sum(true_deact));
        ROPE_only(1).null_incorr(i,j) = 100*ROPE_only_incorr_null/(total - sum(true_null)); 
        HDI_ROPE(1).pos_incorr(i,j)   = 100*HDI_ROPE_incorr_pos/(total - sum(true_act));
        HDI_ROPE(1).neg_incorr(i,j)   = 100*HDI_ROPE_incorr_neg/(total - sum(true_deact)); 
        HDI_ROPE(1).null_incorr(i,j)  = 100*HDI_ROPE_incorr_null/(total - sum(true_null));
        
        clearvars -except ROPE_only HDI_ROPE NHST prior_SD sim_path Nstep rep con_files_folder binary_mask i 
    end
end