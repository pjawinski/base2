% =============================================================================================
% === Permutation-based analysis to test for stronger effects than expected under the null  ===
% =============================================================================================
% /opt/matlab/bin/matlab -nodesktop -nodisplay 

% set working directory
cd '/slow/projects/base2'

% add path with custom functions
addpath code/functions

% load data
clear all
master = importdata('data/03_brainage_w_phenotypes.txt');
master.varnames = master.textdata(1,2:size(master.textdata,2))';
master.IID = master.textdata(2:size(master.textdata,1),1);
master = rmfield(master,'textdata')

% get variables of interest
index = find(startsWith(master.varnames,{'brainage_gap_gm_stack','brainage_gap_wm_stack','brainage_gap_gwm_stack'}));
brainAGE = master.data(:,index');
master.varnames(index')

covs = struct();
index = find(startsWith(master.varnames,{'sex','age','age2', 'TIV'}));
covs.data = master.data(:,index');
covs.varnames = master.varnames(index',1);
covs.varnames

vars = struct();
dep = {'Educ_final',
  'hnetto',
  'MMSE_Summe',
  'GDS_Summe',
  'CESD_Summe',
  'Rauchen_aktuell_inverted',
  'Alkohol_haufigkeit',
  'Alkohol_Menge',
  'Alkohol_6Glaser',
  'CH_Diabetes',
  'HOMAIR',
  'HbA1c',
  'BZP1',
  'BZP2',
  'BMI',
  'RRdi',
  'RRsy',
  'finalMetLscore',
  'GammaGTGGTUL',
  'HarnsaeuremgdL',
  'TNF1',
  'DS2_corr',
  'EM_final',
  'WM_final',
  'Gf_final',
  'futi_mean',
  'cfc_mean'}

index = NaN(size(dep,1),1);
for i = 1:size(dep,1)
    index(i) = find(startsWith(master.varnames,dep(i)));
end
 
vars.data = master.data(:,index');
vars.varnames = master.varnames(index',1);
vars.varnames

%% hypotheses: set expected effect directions
hypo = [...
            0   %'Educ_final'
            0   %'hnetto'
            0   %'MMSE_Summe'
            1   %'GDS_Summe'
            1   %'CESD_Summe'
            1   %'Rauchen_aktuell_inverted'
            1   %'Alkohol_haufigkeit'
            1   %'Alkohol_Menge'
            1   %'Alkohol_6Glaser'
            1   %'CH_Diabetes'
            1   %'HOMAIR'
            1   %'HbA1c'
            1   %'BZP1'
            1   %'BZP2'
            1   %'BMI'
            1   %'RRdi'
            1   %'RRsy'
            1   %'finalMetLscore'
            1   %'GammaGTGGTUL'
            1   %'HarnsäuremgdL'
            1   %'TNF1'
            0   %'DS2_corr'
            0   %'EM_final'
            0   %'WM_final'
            0   %'Gf_final'
            0   %'futi_mean
            0]; %'cfc_mean'
  
% calculate correlations
[rho_obs_partial_original, ~] = partialcorr(brainAGE, vars.data, covs.data, 'rows', 'pairwise', 'type', 'Pearson');
  
% calculate the proportion of tests with hypothesis-consistent effect
% directions - 70.37%, 85.19%, 85.19%!
sum([rho_obs_partial_original' > 0] == hypo)/size(hypo,1)
hits_obs = sum([rho_obs_partial_original' > 0] == hypo)

% invert variables with expected negative directions to facilitate
% one-tailed tests
vars.data_transformed = zeros(size(vars.data));
invert = hypo;
invert(hypo==0) = -1;
vars.data_inv = vars.data.*invert';

% run one-tailed tests
[rho_obs_partial, p_obs_partial] = partialcorr(brainAGE, vars.data_inv, covs.data, 'rows', 'pairwise', 'tail', 'right');

% get pairwise n 
n = NaN(size(vars.data_inv,2),1);
for i = 1:size(vars.data_inv,2)
    n(i,1) = sum(isnan(vars.data_inv(:,i)) == 0);
end

% create output file
table = array2table([n rho_obs_partial_original' p_obs_partial'], 'VariableNames', [{'n'} {'rho_gm'}, {'rho_wm'}, {'rho_gwm'}, {'p_gm'}, {'p_wm'}, {'p_gwm'}], 'RowNames', vars.varnames);
writetable(table, strcat('code/tables/main_corr_MATLAB.txt'), 'Delimiter', '\t', 'WriteRowNames', 1)

% import R results and test whether results are the same - They are!
R = readtable('code/tables/main_corr.txt');
isequal(n, R.n)
isequal(round([table.rho_gm table.rho_wm table.rho_gwm],10),round([R.gm_estimate R.wm_estimate R.gwm_estimate],10))
isequal(round([table.p_gm table.p_wm table.p_gwm],10),round([R.gm_p_value R.wm_p_value R.gwm_p_value],10))

%% Residualize variables before data permutation

% adjust brainage by sex, age, and age2 (due to varying number of observations,
% this needs to be done for each outcome variable seperately).
brainAGE_res = NaN([size(brainAGE) size(vars.data_inv,2)]);
for i = 1:size(brainAGE,2)
    for j = 1:size(vars.data_inv,2)
    logical = isnan(vars.data_inv(:,j))==0;
    model = fitlm(covs.data(logical,:), brainAGE(logical,i));
    brainAGE_res(logical,i,j) = model.Residuals{:,1};
    end
end

% adjust vars by sex, age, and age2 
vars.data_res = NaN(size(vars.data_inv));
for i = 1:size(vars.data_inv,2)
    logical = isnan(vars.data_inv(:,i))==0;
    model = fitlm(covs.data(logical,:), vars.data_inv(logical,i));
    vars.data_res(logical,i) = model.Residuals{:,1};
end

% Test whether correlations of residualized variables correspond to
% partial correlations. - identical with up to 12 decimal places
rho_obs_res = NaN(size(brainAGE,2),size(vars.data_inv,2));
p_obs_res = NaN(size(brainAGE,2),size(vars.data_inv,2));
for i = 1:size(vars.data_inv,2)
    rho_obs_res(:,i) = corr(brainAGE_res(:,:,i), vars.data_res(:,i), 'rows', 'pairwise', 'type', 'Pearson');
    p_obs_res(:,i) = pncovs(rho_obs_res(:,i), sum(isnan(vars.data_inv(:,i))==0),size(covs.data,2), 'right');
end

isequal(round(rho_obs_partial,12),round(rho_obs_res,12)) % should be 1
isequal(round(p_obs_partial,12),round(p_obs_res,12)) % should be 1

% Permutation analysis require a single residualized gm, wm, and gwm
% brainAGE variable. Calculate residuals based on the complete sample and
% correlate results with the outcome-specific brainAGE_res solutions.
% - all correlations above 0.993!
brainAGE_res_multiple = brainAGE_res;
brainAGE_res = NaN(size(brainAGE));
for i = 1:size(brainAGE,2)
    model = fitlm(covs.data, brainAGE(:,i));
    brainAGE_res(:,i) = model.Residuals{:,1};
end

min([corr(brainAGE_res(:,1), reshape(brainAGE_res_multiple(:,1,:), [size(brainAGE_res_multiple,1) size(brainAGE_res_multiple,3)]),'rows', 'pairwise', 'type', 'Pearson')' ...
 corr(brainAGE_res(:,2), reshape(brainAGE_res_multiple(:,2,:), [size(brainAGE_res_multiple,1) size(brainAGE_res_multiple,3)]),'rows', 'pairwise', 'type', 'Pearson')' ...
 corr(brainAGE_res(:,3), reshape(brainAGE_res_multiple(:,3,:), [size(brainAGE_res_multiple,1) size(brainAGE_res_multiple,3)]),'rows', 'pairwise', 'type', 'Pearson')'])

% Measure the absolute discrepencancy between correlations of residualized variables
% (brainAGE residuals based on all subjects) and partial correlations
% - differences in rho and p are negligible
rho_obs_res = NaN(size(brainAGE,2),size(vars.data_inv,2));
p_obs_res = NaN(size(brainAGE,2),size(vars.data_inv,2));
for i = 1:size(vars.data_inv,2)
    rho_obs_res(:,i) = corr(brainAGE_res, vars.data_res(:,i), 'rows', 'pairwise', 'type', 'Pearson');
    p_obs_res(:,i) = pncovs(rho_obs_res(:,i), sum(isnan(vars.data_inv(:,i))==0),size(covs.data,2), 'right');
end

mean(abs(rho_obs_partial' - rho_obs_res')) % mean difference in rho = 0.00007
mean(abs(p_obs_partial' - p_obs_res')) % mean difference in p = 0.00016

max(abs(rho_obs_partial' - rho_obs_res')) % maximum difference in rho = 0.008
max(abs(p_obs_partial' - p_obs_res')) % maximum difference in p = 0.0016

%% Calculate permutation-based p-value for effect directions

% create vector for permutation
npermutations = 1000000;
rng(68882,'twister')

u = zeros(size(brainAGE_res,1), npermutations);
for i = 1:npermutations
    u(:,i) = randperm(size(u,1));
end

% set maximum number of compute threads
maxNumCompThreads(1)

% start parallel pool
myCluster = parcluster('local');
myCluster.NumWorkers = 100;
saveProfile(myCluster);
parpool('local',100);

% Get the number of expected 'hits' (consistent effect directions)
hits_exp = zeros(3, npermutations);
p_exp = zeros(3, size(vars.data_res,2), npermutations);
rho_exp = zeros(3, size(vars.data_res,2), npermutations);

tic
parfor i = 1:npermutations
    [rho_temp, ~ ] = corr(brainAGE_res(u(:,i),:), vars.data_res, 'rows', 'pairwise', 'type', 'Pearson');
    hits_exp(:,i) = sum(rho_temp' > 0);
    p_temp = NaN(size(rho_temp));
    
    for j = 1:size(rho_temp,2)
    p_temp(:,j) = pncovs(rho_temp(:,j), sum(isnan(vars.data_res(:,j))==0),size(covs.data,2), 'right');
    end
    
    rho_exp(:,:,i) = rho_temp;
    p_exp(:,:,i) = p_temp;
end
toc 

% get permutation-based p-values for the observed number of hypothesis-consistent
% effect directions. - 0.1025, 0.0065, 0.0067
% - 0.002 overall
sum(hits_exp' >= hits_obs)/npermutations
sum(sum(hits_exp,1)' >= sum(hits_obs))/(npermutations*3)
sum(hits_obs) % count of associations with hypothesis-consistent effect directions
sum(hits_obs)/size(rho_obs_partial(:),1) % percentage of associations with hypothesis-consistent effect directions

% get permutation-based p-values for the observed number of nominally
% significant results - 0.3469, 0.0010, 0.0027
% - 0.0013 overall
sum(reshape(sum(p_exp < 0.05,2),[size(p_exp,1), size(p_exp,3)])' >= sum(p_obs_partial < 0.05,2)')/npermutations
sum(sum(reshape(p_exp < 0.05,[size(p_exp,1)*size(p_exp,2), size(p_exp,3)]) > 0.05)' >= sum(p_obs_partial(:) < 0.05)')/(npermutations*3)
sum(p_obs_partial(:) < 0.05) % count of nominally significant results
sum(p_obs_partial(:) < 0.05)/size(p_obs_partial(:),1) % percentage of nominally significant results

% get observed p values
rho_obs = rho_obs_partial;
p_obs = p_obs_partial;

% calculate mean observed and expected rho
rho_obs_atanh = atanh(rho_obs);
rho_obs_mean = tanh(mean(rho_obs_atanh,2));
rho_obs_mean_all = tanh(mean(rho_obs_atanh(:)));

rho_exp_atanh = atanh(rho_exp);
rho_exp_mean = tanh(mean(rho_exp_atanh,2));
rho_exp_mean = reshape(rho_exp_mean,[3 npermutations]);
rho_exp_atanh_all = reshape(rho_exp_atanh, [size(rho_exp_atanh,1)*size(rho_exp_atanh,2), size(rho_exp_atanh,3)]);
rho_exp_mean_all = tanh(mean(rho_exp_atanh_all,1));
rho_exp_mean_all = rho_exp_mean_all(:);

% return mean observed rho and respective permutation-based p-values
% 0.0325 (0.0426), 0.0671 (2E-4), 0.0609 (6E-4)
% all: 0.0535 (9E-4)
rho_obs_mean
sum(rho_exp_mean>=rho_obs_mean,2)/npermutations % 1 - (sum(rho_exp_mean>rho_obs_mean,2)/npermutations)

rho_obs_mean_all
sum(rho_exp_mean_all>=rho_obs_mean_all)/npermutations % 1 - (sum(rho_exp_mean_all>rho_obs_mean_all)/npermutations)

% calculate inflation factor 'lambda', i.e. the ratio of the observed median
% chi2-value vs. the expected median chi2-value (chi2inv(0.5,1) = 0.4549)
% Lambda = 2.9352, 5.1468, 4.1591; 3.5179
lambdas = median(chi2inv(1-p_obs_partial',1))/chi2inv(0.5,1)
lambda_all = median(chi2inv(1-p_obs_partial(:)',1))/chi2inv(0.5,1)

    % calculate lambdas from permutations
    p_exp_chi2 = chi2inv(1-p_exp,1);
    p_exp_chi2_gm = reshape(p_exp_chi2(1,:,:), [size(p_exp_chi2,2) size(p_exp_chi2,3)]);
    p_exp_chi2_wm = reshape(p_exp_chi2(2,:,:), [size(p_exp_chi2,2) size(p_exp_chi2,3)]);
    p_exp_chi2_gwm = reshape(p_exp_chi2(3,:,:), [size(p_exp_chi2,2) size(p_exp_chi2,3)]);

    p_exp_chi2_gm_lambdas = median(p_exp_chi2_gm,1)'/chi2inv(0.5,1);
    p_exp_chi2_wm_lambdas = median(p_exp_chi2_wm,1)'/chi2inv(0.5,1);
    p_exp_chi2_gwm_lambdas = median(p_exp_chi2_gwm,1)'/chi2inv(0.5,1);
    p_exp_chi2_lambda_all = median([p_exp_chi2_gm;p_exp_chi2_wm;p_exp_chi2_gwm],1)'/chi2inv(0.5,1);
    
    % Compute the p-values of the observed lambdas as the proportion of
    % permutation-lambdas that exceed the observed lambdas
    % p = 0.0283, 7E-4, 0.0038; 0.0048
    sum(p_exp_chi2_gm_lambdas >= lambdas(1))/npermutations
    sum(p_exp_chi2_wm_lambdas >= lambdas(2))/npermutations
    sum(p_exp_chi2_gwm_lambdas >= lambdas(3))/npermutations
    sum(p_exp_chi2_lambda_all >= lambda_all)/npermutations
   
% how many associations reached p < 0.001 by chance (0.001 = threshold of 
% significance for individual tests after multiple-testing correction)
% - 4.9%
p_exp_all = reshape(p_exp, [size(p_exp,1)*size(p_exp,2) size(p_exp,3)]);
sum(sum(p_exp_all <= 0.001) >= 1)/npermutations

%% Draw one-sided qq-plots for each brainAGE variable

% draw qq plot
ttl{1} = 'Grey matter';
ttl{2} = 'White matter';
ttl{3} = 'Grey and white matter';

qqplot = figure(); hold on;
for j = 1:3
    
    % initiate new subplot
    subplot(1,3,j); hold on;
    set(0,'DefaultTextFontname', 'CMU Serif');
    set(0,'defaulttextinterpreter','latex')
    set(gca,'TickLabelInterpreter', 'latex');

    % sort p-values derived from each permutation and convert to Z-values
    RND_Vector_p = reshape(sort(p_exp(j,:,:), 'descend'),[size(p_exp,2) size(p_exp,3)]);
    RND_Vector_Z = atanh(RND_Vector_p); % return p values again by: RND_Vector_Z_p = 2*normcdf(RND_Vector_Z);

    % get mean expected p-values and confidence intervals from permutation-based
    % p-values at rank 1:ntests
    Mean = zeros(size(RND_Vector_Z,1),1);
    Konf95 = zeros(size(RND_Vector_p,1),1);
    Konf05 = zeros(size(RND_Vector_p,1),1);

    for i=1:size(RND_Vector_p,1)
        Mean(i,1) = -log10(tanh(mean(RND_Vector_Z(i,:))));    
        Konf95(i,1) = -log10(min(maxk(RND_Vector_p(i,:),size(RND_Vector_p,2)/20)));
        Konf05(i,1) = -log10(max(mink(RND_Vector_p(i,:),size(RND_Vector_p,2)/20)));
    end

    % plot title
    title(ttl{j},'fontsize',10,'Interpreter','latex');

    % draw confidence interval
    x = Mean;
    x2 = [x;flipud(x)];
    inBetween = [Konf05;flipud(Konf95)];
    a = fill(x2, inBetween, [211,211,211]/255);
    set(a,'EdgeColor','none','facealpha',0.5);
    
    % plot observed p-values against expected p-values (log-scale)
    plot(Mean,-log10(sort(p_obs(j,:),'descend')),'o','Color',[0, 0.4470, 0.7410],'Markersize',2,'LineWidth',0.25)

    % plot diagonal line of expected p-values
    plot(Mean,Mean,'k-');

    % set plot styles
    ax = gca;
    ax.XAxis.FontSize = 8;
    ax.YAxis.FontSize = 8;
    ax.TickLabelInterpreter='latex';

    ylim([0 4]); xlim([0 2]);
    xlabel(sprintf('expected p-value ($-log_{10}$ scale)'),'fontsize',10,'Interpreter','latex');
    ylabel(sprintf('observed p-value ($-log_{10}$ scale)'),'fontsize',10,'Interpreter','latex');

    xticks(0:1:3);
    xticklabels(0:1:2);
    ytickformat('%.0f');
    yticks(0:1:4);
    yticklabels(0:1:4);
end

% finish plot
qqplot.Renderer = 'Painter';
set(qqplot,'PaperUnits', 'centimeters');
set(qqplot,'PaperPosition', [0 19 20 8]);
print('code/figures/qqplot.pdf','-dpdf');

%% Draw qq-plot for all tested associations (3 brainAGE x 27 outcome variables)

% sort p-values derived from each permutation and convert to Z-values
RND_Vector_p = sort(reshape(p_exp,[size(p_exp,1)*size(p_exp,2), size(p_exp,3)]), 'descend');
RND_Vector_Z = atanh(RND_Vector_p);

% get mean expected p-values and confidence intervals from permutation-based
% p-values at rank 1:ntests
Mean = zeros(size(RND_Vector_p,1),1);
Konf95 = zeros(size(RND_Vector_p,1),1);
Konf05 = zeros(size(RND_Vector_p,1),1);

for i=1:size(RND_Vector_p,1)
    Mean(i,1) = -log10(tanh(mean(RND_Vector_Z(i,:))));    
    Konf95(i,1) = -log10(min(maxk(RND_Vector_p(i,:),size(RND_Vector_p,2)/20)));
    Konf05(i,1) = -log10(max(mink(RND_Vector_p(i,:),size(RND_Vector_p,2)/20)));
end

% initiate new plot
qqplot = figure(); hold on;
set(0,'DefaultTextFontname', 'CMU Serif');
set(0,'defaulttextinterpreter','latex')
set(gca,'TickLabelInterpreter', 'latex');

    % draw confidence interval
    x = Mean;
    x2 = [x;flipud(x)];
    inBetween = [Konf05;flipud(Konf95)];
    a = fill(x2, inBetween, [211,211,211]/255);
    set(a,'EdgeColor','none','facealpha',0.5);

    % plot observed p-values against expected p-values (log-scale)
    plot(Mean,-log10(sort(p_obs(:),'descend')),'o','Color',[0, 0.4470, 0.7410],'Markersize',2,'LineWidth',0.25)

    % plot diagonal line of expected p-values
    plot(Mean,Mean,'k-');

    % set plot styles
    ax = gca;
    ax.XAxis.FontSize = 8;
    ax.YAxis.FontSize = 8;
    ax.TickLabelInterpreter='latex';

    ylim([0 4]); xlim([0 2]);
    xlabel(sprintf('expected p-value ($-log_{10}$ scale)'),'fontsize',10,'Interpreter','latex');
    ylabel(sprintf('observed p-value ($-log_{10}$ scale)'),'fontsize',10,'Interpreter','latex');

    xticks(0:1:3);
    xticklabels(0:1:2);
    ytickformat('%.0f');
    yticks(0:1:4);
    yticklabels(0:1:4);

% finish plot
qqplot.Renderer = 'Painter';
set(qqplot,'PaperUnits', 'centimeters');
set(qqplot,'PaperPosition', [2 18 5 8]);
print('code/figures/qqplot_all.pdf','-dpdf');

%% delete parcluster
p = gcp;
delete(p);
myCluster = parcluster('local');
delete(myCluster.Jobs)
