function generate_con_files(sim_path,iXYZ,act,deact,background,mask_file,PSC,SD,sk,ku,PSC_bg,N,XYZ)

tic
% activaions, deactivations and practically non-significant background
signal = PSC*act - PSC*deact + PSC_bg*background.*(2*randi([0 1],1,length(background))-1);
noise = pearsrnd(0,SD,sk,ku,[N, length(signal)]); 

hdr = spm_vol(mask_file);
mask = spm_read_vols(hdr);
signal_img = mask;
signal_img(iXYZ) = signal;

status = exist([sim_path filesep 'sim_con_files' filesep 'PSC_(' num2str(PSC) ')_SD_(' num2str(SD) ')_Sk_(' num2str(sk) ')_Ku_(' num2str(ku) ')']);
if status ~= 7
   mkdir([sim_path filesep 'sim_con_files' filesep 'PSC_(' num2str(PSC) ')_SD_(' num2str(SD) ')_Sk_(' num2str(sk) ')_Ku_(' num2str(ku) ')']); 
end

for i=1:N
    hdr.fname = [sim_path filesep 'sim_con_files' filesep 'PSC_(' num2str(PSC) ')_SD_(' num2str(SD) ')_Sk_(' num2str(sk) ')_Ku_(' num2str(ku) ')' filesep 'con_' num2str(i,'%04d') '.nii'];
    noise_img    = mask;
    noise_img(iXYZ) = noise(i,:);
    con_img = signal_img + noise_img;  
    spm_write_vol(hdr,con_img);
    clear noise_img con_img scon_img
end

time = toc;
fprintf(['Generate :: PSC_(' num2str(PSC) ')_SD_(' num2str(SD) ')_Sk_(' num2str(sk) ')_Ku_(' num2str(ku) ') :: Done in ' num2str(time)  filesep 'n']);

