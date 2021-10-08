% ===============================================================================
% === Plot brain-predicted vs. chronological age with matched ukb individuals ===
% ===============================================================================

% set working directory
cd '/Users/philippe/Desktop/base2'

% import data
data = importdata('code/derivatives/04_matched_datasets.txt');
data.varnames = data.textdata(1,2:10);
data.IID = data.textdata(2:end,1);
data = rmfield(data,'textdata');

% load ukb and base2 data
primary = data;
select = primary.data(:,1) == 1;
primary.data = primary.data(select,:);
primary.IID = primary.IID(select,:);

secondary = data;
select = secondary.data(:,1) == 2;
secondary.data = secondary.data(select,:);
secondary.IID = secondary.IID(select,:);

% calculate correlations and mean absolute errors
primary.age = primary.data(:,3);
corr(primary.age, primary.data(:,5:7))
mean(abs(primary.age - primary.data(:,5:7)))

secondary.age = secondary.data(:,3);
corr(secondary.age, secondary.data(:,5:7))
mean(abs(secondary.age - secondary.data(:,5:7)))

% make directory for figures
system('mkdir -p code/figures/');

% calculate correlations and mean absolute errors
ttl = {'Grey matter','White matter','Grey and white matter'};
txt = cell(1);
txt2 = cell(1);
for i = 1:3
    txt{i} = compose(strcat({'MAE = '}, num2str(round(mean(abs(primary.age - primary.data(:,i+4))),2),'%1.2f'), {' yrs\n rho = '}, num2str(round(corr(primary.age, primary.data(:,i+4)),3),'%1.3f')));
    txt2{i} = compose(strcat({'MAE = '}, num2str(round(mean(abs(secondary.age - secondary.data(:,i+4))),2),'%1.2f'), {' yrs\n rho = '}, num2str(round(corr(secondary.age, secondary.data(:,i+4)),3),'%1.3f')));
end

% plot brain-predicted vs. chronological age (stacked models)
brainage_figure(1) = figure(); hold on;
for i = 1:3
subplot(1,3,i); hold on;
    plot(primary.age, primary.data(:,i+4), '.','Color',[210/255, 210/255, 210/255],'Markersize',2,'LineWidth',0.25)
    plot(secondary.age, secondary.data(:,i+4), '.','Color',[0, 0.4470, 0.7410],'Markersize',2,'LineWidth',0.25)

    set(0,'DefaultTextFontname', 'CMU Serif');
    set(0,'defaulttextinterpreter','latex');
    set(gca,'TickLabelInterpreter', 'latex', 'box','off');

    ax = gca;
    ax.Box = 'off';
    ax.XAxis.FontSize = 7;
    ax.YAxis.FontSize = 7;
    ax.TickLabelInterpreter='latex';

    xlabel(sprintf('chronological age (years)'),'fontsize',9,'Interpreter','latex');
    ylabel(sprintf('brain age (years)'),'fontsize',9,'Interpreter','latex');
    
    ylim([42 85]); xlim([42 85]);
    %ylim([53 85]); xlim([53 85]);
    xticks(45:10:85);
    xticklabels(45:10:85);
    yticks(45:10:85);
    yticklabels(45:10:85);

    % axis label position
    xh = get(gca,'xlabel');
    get(xh,'position');
    set(xh,'position', [63.5 37.5 -1]); 
    
    % prevent x axis cut off
    pos = get(gca, 'Position');
    set(gca, 'Position', [pos(1) pos(2)+0.05 pos(3) pos(4)-0.05]);

    text(66,48,txt{i},'Fontsize', 7, 'HorizontalAlignment', 'left');
    text(45,78,txt2{i},'Fontsize', 7, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);
    title(ttl{i}, 'Fontsize', 8);
    
end

set(brainage_figure(1),'PaperUnits', 'centimeters');
set(brainage_figure(1),'PaperPosition', [0.5 20 20 6.5]);
print('code/figures/accuracy_matched.pdf','-dpdf')
print('code/figures/accuracy_matched.png','-dpng', '-r300')
