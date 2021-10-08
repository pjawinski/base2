% =======================
% === Read CAT12 data ===
% =======================
% /opt/matlab/bin/matlab -nodesktop -nodisplay 

% set working directory
cd '/slow/projects/base2/'

% add paths
addpath /slow/projects/base2/code/functions/
addpath /fast/software/matlab/spm12/
addpath /fast/software/matlab/RVM/

% set maximum number of threads used for computations
maxNumCompThreads(28);

% find _cat12.zip files
files = dir('data/T1w/*/*_cat12.zip');

% get individual identifier from files.name
IID = {};
zeros(numel(files),1);
for i = 1:numel(files)
    temp = split(files(i).name, '_');
    IID{i,1} = temp{2};
end

% unzip rsmwp[1/2]T1_orig_defaced.nii, read nifti, and create grey and white matter matrix
files_constant = parallel.pool.Constant(files);
gm = zeros(numel(files), 16128);
wm = zeros(numel(files), 16128);

tic
parfor i = 1:numel(files)
    files_c = files_constant.Value; 
    system(char(strcat({'unzip -joqq '}, files(i).folder, '/', files_c(i).name, {' cat12/mri/rsmwp1*.nii -d '}, files_c(i).folder, '/')));
    rsmwp1 = dir(strcat(files(i).folder, '/rsmwp1*.nii'))
    nifti = read_nifti(char(strcat(files_c(i).folder, '/', rsmwp1.name)));
    gm(i,:) = nifti(:)';
    system(char(strcat({'rm -f '}, files_c(i).folder, '/', rsmwp1.name))); 
    rsmwp1 = [];
    nifti = [];
    
    system(char(strcat({'unzip -joqq '}, files_c(i).folder, '/', files_c(i).name, {' cat12/mri/rsmwp2*.nii -d '}, files_c(i).folder, '/')));
    rsmwp2 = dir(strcat(files(i).folder, '/rsmwp2*.nii'))
    nifti = read_nifti(char(strcat(files_c(i).folder, '/', rsmwp2.name)));
    wm(i,:) = nifti(:)';
    system(char(strcat({'rm -f '}, files_c(i).folder, '/', rsmwp2.name))); 
    rsmwp2 = [];
    nifti = [];
end
toc

% test if all files have been imported. - Yep.
sum(sum(gm,2) == 0) % should be 0
sum(sum(wm,2) == 0) % should be 0

% clean work space
clearvars -except gm files wm parent IID

% get features from UKB dataset
load('code/models/RVM_gm.mat')
load('code/models/RVM_wm.mat')
gm = gm(:,RVM_gm.train_logical);
wm = wm(:,RVM_wm.train_logical);

%%
% ===================
% === Predict Age ===
% ===================

% -------
% xgboost 
% -------

% gm: transform to pca scores
Test_Samples = gm;
x = load('code/models/xgb_gm_pca_training.mat');
Test_Samples_centered = Test_Samples-repmat(x.Train_Samples_means,size(Test_Samples,1),1);
Test_Samples_pca_score = Test_Samples_centered/transpose(x.Train_Samples_pca_coeff);
dlmwrite('code/models/xgb_gm_Test_Samples_pca.txt', Test_Samples_pca_score, 'delimiter', '\t', 'precision', 16)
 
% wm: transform to pca scores
Test_Samples = wm;
x = load('code/models/xgb_wm_pca_training.mat');
Test_Samples_centered = Test_Samples-repmat(x.Train_Samples_means,size(Test_Samples,1),1);
Test_Samples_pca_score = Test_Samples_centered/transpose(x.Train_Samples_pca_coeff);
dlmwrite('code/models/xgb_wm_Test_Samples_pca.txt', Test_Samples_pca_score, 'delimiter', '\t', 'precision', 16)
 
% predict age in R, import age estimates, clean up
system('code/models/xgboost_pred.R');
pred_xgb = importdata('code/models/xgb_pred.txt');
delete('code/models/*.txt')

% -------
% rvm
% -------

% get rvm_gm estimates
Test_Samples = gm;
RVM_gm.rv_index = 1:size(RVM_gm.rv_mu,1);
Test_Samples_centered = Test_Samples-repmat(RVM_gm.train_means,size(Test_Samples,1),1);
Test_Samples_pca_score = Test_Samples_centered/RVM_gm.train_pca_coeff';
[gm_y_mu, ~] = rvm_test(RVM_gm,Test_Samples_pca_score);

% get rvm_wm estimates
Test_Samples = wm;
RVM_wm.rv_index = 1:size(RVM_wm.rv_mu,1);
Test_Samples_centered = Test_Samples-repmat(RVM_wm.train_means,size(Test_Samples,1),1);
Test_Samples_pca_score = Test_Samples_centered/RVM_wm.train_pca_coeff';
[wm_y_mu, ~] = rvm_test(RVM_wm,Test_Samples_pca_score);
pred_rvm = [gm_y_mu wm_y_mu];

% ------------------------
% merge into one structure
% ------------------------

% create structure
brainage = struct();
brainage.data = NaN(size(gm,1),36);
brainage.varnames = {'brainage_gm_xgtree','brainage_gm_xglin', 'brainage_gm_rvm', 'brainage_gm_stack', 'brainage_wm_xgtree', 'brainage_wm_xglin', 'brainage_wm_rvm', 'brainage_wm_stack', 'brainage_gwm_xgtree', 'brainage_gwm_xglin', 'brainage_gwm_rvm', 'brainage_gwm_stack', ...
    'brainage_gap_gm_xgtree','brainage_gap_gm_xglin', 'brainage_gap_gm_rvm', 'brainage_gap_gm_stack', 'brainage_gap_wm_xgtree', 'brainage_gap_wm_xglin', 'brainage_gap_wm_rvm', 'brainage_gap_wm_stack', 'brainage_gap_gwm_xgtree', 'brainage_gap_gwm_xglin', 'brainage_gap_gwm_rvm', 'brainage_gap_gwm_stack', ...
    'brainage_gap_adj_gm_xgtree','brainage_gap_adj_gm_xglin', 'brainage_gap_adj_gm_rvm', 'brainage_gap_adj_gm_stack', 'brainage_gap_adj_wm_xgtree', 'brainage_gap_adj_wm_xglin', 'brainage_gap_adj_wm_rvm', 'brainage_gap_adj_wm_stack', 'brainage_gap_adj_gwm_xgtree', 'brainage_gap_adj_gwm_xglin', 'brainage_gap_adj_gwm_rvm', 'brainage_gap_adj_gwm_stack'};

brainage.data(:,[1 2 5 6]) = pred_xgb;
brainage.data(:,[3 7]) = pred_rvm;

% stack estimates through linear regression
load('code/models/stack_coefficients.mat')
brainage.data(:,4) = stack_coefficients.gm(1) + stack_coefficients.gm(2) * brainage.data(:,1) + stack_coefficients.gm(3) * brainage.data(:,2) + stack_coefficients.gm(4) * brainage.data(:,3);
brainage.data(:,8) = stack_coefficients.wm(1) + stack_coefficients.wm(2) * brainage.data(:,5) + stack_coefficients.wm(3) * brainage.data(:,6) + stack_coefficients.wm(4) * brainage.data(:,7);
brainage.data(:,9) = stack_coefficients.gwm_xgtree(1) + stack_coefficients.gwm_xgtree(2) * brainage.data(:,1) + stack_coefficients.gwm_xgtree(3) * brainage.data(:,5);
brainage.data(:,10) = stack_coefficients.gwm_xglin(1) + stack_coefficients.gwm_xglin(2) * brainage.data(:,2) + stack_coefficients.gwm_xglin(3) * brainage.data(:,6);
brainage.data(:,11) = stack_coefficients.gwm_rvm(1) + stack_coefficients.gwm_rvm(2) * brainage.data(:,3) + stack_coefficients.gwm_rvm(3) * brainage.data(:,7);
brainage.data(:,12) = stack_coefficients.gwm(1) + stack_coefficients.gwm(2) * brainage.data(:,1) + stack_coefficients.gwm(3) * brainage.data(:,2) + stack_coefficients.gwm(4) * brainage.data(:,3) + stack_coefficients.gwm(5) * brainage.data(:,5) + stack_coefficients.gwm(6) * brainage.data(:,6) + stack_coefficients.gwm(7) * brainage.data(:,7);

%% add phenotypic data

% import demographics
demographics = importdata('code/derivatives/01_phenotypes_master.txt');
demographics.varnames = demographics.textdata(1,2:size(demographics.textdata,2))';
demographics.IID = demographics.textdata(2:size(demographics.textdata,1),1);
demographics = rmfield(demographics,'textdata');

% test whether gm/wm IID is aligned with demographics % FALSE
isequal(IID, demographics.IID)

% align data
[~, IA, IB] = intersect(IID, demographics.IID);
isequal(IID(IA),demographics.IID(IB))

IID = IID(IA);
brainage.data = brainage.data(IA,:);
gm = gm(IA,:);
wm = wm(IA,:);
files = files(IA);
demographics.IID = demographics.IID(IB);
demographics.data = demographics.data(IB,:);

% test whether gm/wm IID is aligned with demographics - TRUE
isequal(IID,demographics.IID)

%% get quality ratings

% unzip and get ratings
qc_ratings = zeros(length(files),9);

tic
for i = 1:length(files)
    %files_c = files_constant.Value;  
    system(char(strcat({'unzip -joqq '}, files(i).folder, '/', files(i).name, {' cat12/report/cat*.mat -d '}, files(i).folder, '/')));
    cat_file = dir(strcat(files(i).folder, '/cat*.mat'));
    load(strcat(files(i).folder, '/', cat_file.name)); 
    qc_ratings(i,:) = [S.subjectmeasures.vol_TIV, S.subjectmeasures.vol_abs_CGW(:,1:4) S.qualityratings.res_RMS S.qualityratings.NCR S.qualityratings.ICR S.qualityratings.IQR];
    clearvars S
    system(char(strcat({'rm -f '}, files(i).folder, '/', cat_file.name)));
    fprintf(strcat(sprintf('%0.2f',i/length(files)*100), ' percent done.\n'));
end
toc

% have all files been successfully opened?
sum(sum(qc_ratings, 2)==0) % should be 0

% merge variable
qc = struct();
qc.files = files;
qc.ratings = qc_ratings;
qc.ratings_varnames = {'TIV', 'CSF', 'GM', 'WM', 'RES', 'RMS','NCR','ICR', 'IQR'};
qc.IID = IID;

%% select individuals

% remove individuals with no sex/age info, age < 40, or poor IQR (> 3)
selection = sum([isnan(demographics.data(:,1:2)) (demographics.data(:,2) < 40) (qc.ratings(:,9)>3)],2) < 1;
index = find(selection==0)'

demographics.data(index,:) = [];
demographics.IID(index,:) = [];
IID(index,:) = [];
brainage.data(index,:) = [];
gm(index,:) = [];
wm(index,:) = [];
files(index,:) = [];
qc.IID(index,:) = [];
qc.files(index,:) = [];
qc.ratings(index,:) = [];

%% calculate brainage gap

% calculate gap
sex = demographics.data(:,1);
age = demographics.data(:,2);

for i = 1:12
    brainage.data(:,i+12) = brainage.data(:,i) - age;
end

% adjust brainage gap by sex, age, age2, TIV
age2 = demographics.data(:,2).^2;
TIV = qc.ratings(:,1);

for i = 1:12
    model = fitlm([sex age age2 TIV], brainage.data(:,i+12));
    brainage.data(:,i+24) = model.Residuals.Raw;
end

% compare mae and correlations of brainAGE with chronological age
[corr(age, brainage.data(:,1:12))' ...
mean(abs(brainage.data(:,13:24)))']

% save variables
gm_logical = RVM_gm.train_logical;
wm_logical = RVM_wm.train_logical;
save('code/derivatives/03_brainage.mat', 'IID', 'brainage', 'qc', 'gm', 'wm', 'files', 'gm_logical', 'wm_logical', 'demographics')

%% make text files with brain age and phenotype variables

% convert to table
rowNames = IID;
colNames = ['sex', 'age', 'age2', qc.ratings_varnames, brainage.varnames];
brainage_table = array2table([demographics.data(:,1:2), demographics.data(:,2).^2, qc.ratings, brainage.data],'RowNames',rowNames,'VariableNames',colNames);

% write table
writetable(brainage_table, 'code/derivatives/03_brainage.txt', 'Delimiter', '\t', 'WriteRowNames', 1)

% create table with phenotypes
rowNames = IID;
colNames = ['sex', 'age', 'age2', demographics.varnames(3:end)', qc.ratings_varnames, brainage.varnames];
brainage_table = array2table([demographics.data(:,1:2), demographics.data(:,2).^2, demographics.data(:,3:end), qc.ratings, brainage.data],'RowNames',rowNames,'VariableNames',colNames);

% write table
writetable(brainage_table, 'code/derivatives/03_brainage_w_phenotypes.txt', 'Delimiter', '\t', 'WriteRowNames', 1)
