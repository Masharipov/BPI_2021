function bayesian_estimation(sim_path,N,Nstep,rep,con_files_folder,binary_mask)

tic

%% Random sampling
for i = 1:length(Nstep)
    for j = 1:rep
        random_sub(i,j).perm = randperm(N,Nstep(i));
    end
end

%% Estimation
f = waitbar(0,['Estimation']);
for i = 1:length(Nstep) 
    for j = 1:rep
        matlabbatch{1}.spm.stats.factorial_design.dir = {[sim_path filesep 'group_stat' filesep con_files_folder filesep 'Sample_' num2str(Nstep(i),'%03d') '_Rep_' num2str(j,'%03d')]};
        for k = 1:Nstep(i)
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans(k,1) = ...
            {[sim_path filesep 'sim_con_files' filesep con_files_folder filesep 'con_' num2str(random_sub(i,j).perm(k),'%04d') '.nii,1']};
        end
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {binary_mask};
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
        clear matlabbatch
    end
    waitbar(i/length(Nstep),f,'Estimation')
end

delete(f)
save([sim_path filesep 'group_stat' filesep con_files_folder filesep 'random_list_of_subjects.mat'],'random_sub');

time = toc;
fprintf(['Estimate :: ' con_files_folder ' :: Done in ' num2str(time) filesep 'n']);