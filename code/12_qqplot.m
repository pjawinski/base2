% ====================
% == Draw qq-plots ===
% ====================

% set working directory
cd '/slow/projects/base2'

% load data
load('code/derivatives/11_permutations.mat', 'p_exp', 'p_obs');

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
print('code/figures/qqplot.png','-dpng', '-r300')

%% Draw merged qq-plot for all tested associations (3 brainAGE x 27 outcome variables)

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
print('code/figures/qqplot_all.png','-dpng', '-r300')
