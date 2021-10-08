% ===========================
% === CAT12 Preprocessing ===
% ===========================
% /opt/matlab/bin/matlab -nodesktop -nodisplay 

% set working directory
cd '/slow/projects/base2/'

% set mail account to submit error messages 
E_mail = 'account@host.de';
SMTP_Server = 'mailhost.host.de';
SMTP_Username = 'account@host.de';
SMTP_Password = 'pw';
E_mail_target = 'account@host.de';

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
setpref('Internet','E_mail',E_mail);
setpref('Internet','SMTP_Server',SMTP_Server);
setpref('Internet','SMTP_Username',SMTP_Username);
setpref('Internet','SMTP_Password',SMTP_Password);

% set Matlab toolbox paths
addpath /slow/projects/base2/code/functions/
addpath /fast/software/matlab/spm12/

% set temporary folder
path = strcat(pwd,'/');
temppath = strcat(path,'temp/');
system(char(strcat({'mkdir -p '}, temppath)));

% get files that require cat12 preprocessing
scans = dir('data/T1w/**/*.nii')
clearvars scans2remove
k = 0
for i = 1:size(scans,1)
    if exist(strcat(scans(i).folder,'/',scans(i).name(1:(end-4)),'_cat12.zip')) == 2
        k = k + 1;
        scans2remove(k,1) = i;
    end
end
if exist('scans2remove', 'var') == 1
    scans(scans2remove) = []
end

% get subs and folder
subs = extractBefore({scans.name}', '.');
subs_folder = {scans.folder}';

% set up Matlab job (cat12.5 build 1364)
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 0;
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP = 1070;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.darteltpm = {'/fast/software/matlab/spm12/toolbox/cat12/templates_1.50mm/Template_1_IXI555_MNI152.nii'};
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shootingtpm = {'/fast/software/matlab/spm12/toolbox/cat12/templates_1.50mm/Template_0_IXI555_MNI152_GS.nii'};
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regstr = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.fixed = [1 0.1];
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 2;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 2;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobian.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [0 0];
matlabbatch{2}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';

% start parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = 100;
saveProfile(myCluster);
parpool('local',100);

% submit constant variables to workers
matlabbatch_constant = parallel.pool.Constant(matlabbatch);
subs_constant = parallel.pool.Constant(subs);
subs_folder_constant = parallel.pool.Constant(subs_folder);
temppath_constant = parallel.pool.Constant(temppath);

% initiate spm
spm('Defaults','fMRI');
spm_jobman('initcfg');

parfor idx = 1:numel(subs)
    
    % get vars
    temp = temppath_constant.Value; % temp = temppath;
    sub = string(subs_constant.Value(idx)); % sub = string(subs(idx));
    sub_folder = string(subs_folder_constant.Value(idx)); % sub_folder = subs_folder(idx);
    matlabbatch_loop = matlabbatch_constant.Value; % matlabbatch_loop = matlabbatch;

    try
        % copy file to temporary folder and gunzip .nii.gz for spm
        system(strcat({'mkdir -p '}, temp, sub, {'; cp -R '}, sub_folder, {'/. '}, temp, sub));
        % system(strcat({'gunzip -fk '}, temp, sub, '/', sub, '.nii.gz'));    

        if exist(strcat(temp, sub, '/', sub, '.nii'), 'file') == 2
            system(strcat({'mkdir -p '}, temp, sub, '/cat12'));
            system(strcat({'mv '}, temp, sub, '/', sub, {'.nii '}, temp, sub, '/cat12'));

            matlabbatch_loop{1}.spm.tools.cat.estwrite.data = {char(strcat(temp, sub, '/cat12/', sub, '.nii'))};
            matlabbatch_loop{2}.spm.spatial.smooth.data =  {char(strcat(temp, sub, '/cat12/mri/mwp1', sub, '.nii'))
                                                        char(strcat(temp, sub, '/cat12/mri/mwp2', sub, '.nii'))
                                                        char(strcat(temp, sub, '/cat12/mri/wm', sub, '.nii'))};
            spm_jobman('run',matlabbatch_loop);

            resize_img(char(strcat(temp, sub,'/cat12/mri/smwp1', sub, '.nii')), [8 8 8], nan(2,3))
            resize_img(char(strcat(temp, sub,'/cat12/mri/smwp2', sub, '.nii')), [8 8 8], nan(2,3))
            resize_img(char(strcat(temp, sub,'/cat12/mri/swm', sub, '.nii')), [8 8 8], nan(2,3))

            system(strcat({'rm '}, temp, sub, '/cat12/', sub, '.nii'));
            system(strcat({'cd '}, temp, sub, {'/; zip -r '}, sub, '_cat12.zip cat12/ -x "*.DS_Store"'));
            system(strcat({'cp '}, temp, sub, '/', sub, {'_cat12.zip '}, sub_folder, '/'));
        end
        system(strcat({'rm -R '}, temp, sub));    
        
    catch ME
    sendmail(E_mail_target,'Message from Cluster', string(strcat(sub,{': '}, ME.message))); 
    end
end 

sendmail(E_mail_target,'Message from Cluster','Loop finished.');

% shut down cluster and remove all jobs created with profile local
p = gcp;
delete(p);
myCluster = parcluster('local');
delete(myCluster.Jobs)
