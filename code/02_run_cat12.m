% ===========================
% === CAT12 Preprocessing ===
% ===========================
% /opt/matlab/bin/matlab -nodesktop -nodisplay 

% set working directory
cd /slow/projects/base2/data/T1w

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

path = strcat(pwd,'/');
temppath = strcat(path,'temp/');
system(char(strcat({'mkdir -p '}, temppath)));
subjects = importdata(strcat(path,'subs_to_preprocess.txt'));

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

% set maximum number of compute threads
maxNumCompThreads(1)

% start parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = 100;
saveProfile(myCluster);
parpool('local',100);

% submit constant variables to workers
matlabbatch_constant = parallel.pool.Constant(matlabbatch);
subjects_constant = parallel.pool.Constant(subjects);
dir_constant = parallel.pool.Constant(path);
tempdir_constant = parallel.pool.Constant(temppath);

% initiate spm
spm('Defaults','fMRI');
spm_jobman('initcfg');

% loop
parfor idx = 1:size(subjects,1)
    
    dir = dir_constant.Value; % dir = path;
    tempdir = tempdir_constant.Value; % tempdir = temppath;
    element = char(string(subjects_constant.Value(idx))); % element = char(string(subjects(idx)));
    matlabbatch_loop = matlabbatch_constant.Value; % matlabbatch_loop = matlabbatch;
 
    folder = element(4:10);
    
    try
        system(strcat({'cp -R '}, dir, string(folder), {' '}, tempdir));

        if exist(strcat(tempdir, string(folder), '/', string(element)),'file') == 2
            system(strcat({'mkdir '}, tempdir, string(folder), '/cat12'));
            system(strcat({'mv '}, tempdir, string(folder), '/', string(element), {' '}, tempdir, string(folder), '/cat12/'));

            matlabbatch_loop{1}.spm.tools.cat.estwrite.data = {char(strcat(tempdir, string(folder), '/cat12/', string(element)))};
            matlabbatch_loop{2}.spm.spatial.smooth.data = {char(strcat(tempdir, string(folder),'/cat12/mri/mwp1', string(element)))
                                                           char(strcat(tempdir, string(folder),'/cat12/mri/mwp2', string(element)))
                                                           char(strcat(tempdir, string(folder),'/cat12/mri/wm', string(element)))};
            
            cd(strcat(strcat(tempdir, string(folder), '/cat12/')))
            spm_jobman('run',matlabbatch_loop);

            resize_img(char(strcat(tempdir, string(folder),'/cat12/mri/smwp1', string(element))), [8 8 8],nan(2,3))
            resize_img(char(strcat(tempdir, string(folder),'/cat12/mri/smwp2', string(element))), [8 8 8],nan(2,3))
            resize_img(char(strcat(tempdir, string(folder),'/cat12/mri/swm', string(element))), [8 8 8],nan(2,3))

            system(strcat({'rm '}, tempdir, string(folder), '/cat12/', string(element)));
            system(strcat({'cd '}, tempdir, string(folder), {'/; zip -r '}, string(folder), '_cat12.zip cat12/ -x "*.DS_Store"'));

            system(strcat({'cp '}, tempdir, string(folder), '/', string(folder), {'_cat12.zip '}, dir, string(folder), '/'));
        else
            sendmail(E_mail_target,'Message from Cluster', string(strcat(string(folder),{': '}, string(element), {' does not exist'})));
        end

        system(strcat({'rm -R '}, tempdir, string(folder)));
        system(strcat('sed -i -e "s%', string(element), '%%g"', {' '}, strcat(dir,'subs_to_preprocess.txt'))); % remove subj from to do list
        system(string(strcat('sed -i -e "/^[[:space:]]*$/d"', {' '}, strcat(dir,'subs_to_preprocess.txt')))); % remove nasty white space
        
    catch ME
    sendmail(E_mail_target,'Message from Cluster', string(strcat(string(folder),{': '}, ME.message))); 
    end 
end

sendmail(E_mail_target,'Message from Cluster','Loop finished.');

% shut down cluster and remove all jobs created with profile local
p = gcp;
delete(p);
myCluster = parcluster('local');
delete(myCluster.Jobs)
