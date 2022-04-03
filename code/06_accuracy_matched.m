% ===============================================================================
% === Plot brain-predicted vs. chronological age with matched ukb individuals ===
% ===============================================================================

% set working directory
cd '/Users/philippe/desktop/projects/base2'

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

% calculate correlations, mean absolute errors, and weighted mean absoluted
% error
primary.age = primary.data(:,3);
primary.agerange = (max(primary.age)-min(primary.age));
primary.rho = corr(primary.age, primary.data(:,5:7));
primary.mae = mean(abs(primary.age - primary.data(:,5:7)));
primary.wmae = mean(abs(primary.age - primary.data(:,5:7)))/primary.agerange;

secondary.age = secondary.data(:,3);
secondary.agerange = (max(secondary.age)-min(secondary.age));
secondary.rho = corr(secondary.age, secondary.data(:,5:7));
secondary.mae = mean(abs(secondary.age - secondary.data(:,5:7)));
secondary.wmae = mean(abs(secondary.age - secondary.data(:,5:7)))/secondary.agerange;

% make directory for figures
system('mkdir -p code/figures/');

% calculate correlations and mean absolute errors
ttl = {'Grey matter','White matter','Grey and white matter'};
primary.txt = cell(1);
secondary.txt = cell(1);
for i = 1:3
    primary.txt{i} = compose(strcat({'MAE = '}, num2str(round(primary.mae(i),2),'%1.2f'), {' yrs\nwMAE = '}, num2str(round(primary.wmae(i),2),'%1.2f'), {' yrs\nrho = '},num2str(round(primary.rho(i),3),'%1.3f')));
    secondary.txt{i} = compose(strcat({'MAE = '}, num2str(round(secondary.mae(i),2),'%1.2f'), {' yrs\nwMAE = '}, num2str(round(secondary.wmae(i),2),'%1.2f'), {' yrs\nrho = '},num2str(round(secondary.rho(i),3),'%1.3f')));
end

% plot brain-predicted vs. chronological age (stacked models)
brainage_figure(1) = figure(); hold on;
for i = 1:3
subplot(1,3,i); hold on;

    % scatter plot of brain-predicted vs. chronological age    
    plot(primary.age, primary.data(:,i+4), '.','Color',[210/255, 210/255, 210/255],'Markersize',2,'LineWidth',0.25)
        xlim([min(primary.age) max(primary.age)]); ylim([min(primary.age) max(primary.age)]); 
        mdl = fitlm(primary.age,primary.data(:,i+4));
        prim_refline = refline(mdl.Coefficients{2,1}, mdl.Coefficients{1,1});
        prim_refline.Color = [180/255, 180/255, 180/255];
    
    plot(secondary.age, secondary.data(:,i+4), '.','Color',[0, 0.4470, 0.7410],'Markersize',2,'LineWidth',0.25)
        xlim([min(secondary.age) max(secondary.age)]); ylim([min(secondary.age) max(secondary.age)]); 
        mdl = fitlm(secondary.age,secondary.data(:,i+4));
        sec_refline = refline(mdl.Coefficients{2,1}, mdl.Coefficients{1,1});
        sec_refline.Color = [0, 0.4470, 0.7410];

    % set text interpreter
    set(0,'DefaultTextFontname', 'CMU Serif');
    set(0,'defaulttextinterpreter','latex');
    set(gca,'TickLabelInterpreter', 'latex', 'box','off');
    
    % axis settings
    ylim([42 85]); xlim([42 85]);
    xticks(45:10:85);
    xticklabels(45:10:85);
    yticks(45:10:85);
    yticklabels(45:10:85);
    
    ax = gca;
    ax.Box = 'off';
    ax.XAxis.FontSize = 8;
    ax.YAxis.FontSize = 8;
    ax.TickLabelInterpreter='latex';

    % axis labels and position
    xlabel(sprintf('chronological age (years)'),'fontsize',10,'Interpreter','latex');
    ylabel(sprintf('brain age (years)'),'fontsize',10,'Interpreter','latex');
    xh = get(gca,'xlabel');
    get(xh,'position');
    set(xh,'position', [63.5 37.5 -1]); 
    
    % avoid x axis cut off
    pos = get(gca, 'Position');
    set(gca, 'Position', [pos(1) pos(2)+0.05 pos(3) pos(4)-0.05]);

    % add model accuracy metrics
    text(63.5,79.5,strcat(secondary.txt{i}, '{\hspace{0.4cm}}'), 'Fontsize', 7, 'HorizontalAlignment', 'right', 'Color', [0, 0.4470, 0.7410]);
    text(85,48,strcat(primary.txt{i}, '{\hspace{0.4cm}}'), 'Fontsize', 7, 'HorizontalAlignment', 'right');
    
    % add title
    title(ttl{i}, 'Fontsize', 9);

    % plot x=y identity line
    xy_identity = refline(1,0);
    xy_identity.Color = 'k';
end

set(brainage_figure(1),'PaperUnits', 'centimeters');
set(brainage_figure(1),'PaperSize',[20 6.5]);
set(brainage_figure(1),'PaperPosition', [0 0 20 6.5]);
print('code/figures/accuracy_matched.pdf','-dpdf')
print('code/figures/accuracy_matched.png','-dpng', '-r300')
