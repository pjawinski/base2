% ==================================================
% === Plot brain-predicted vs. chronological age ===
% ==================================================

% set working directory
cd '/users/philippe/desktop/projects/base2'

% load ukb and base2 data
primary = load('data/ukb_brainage.mat');
secondary = load('code/derivatives/03_brainage.mat');

% calculate correlations, mean absolute errors, and weighted mean absoluted
% error
primary.age = primary.covs.data(:,5);
primary.agerange = (max(primary.age)-min(primary.age));
primary.rho = corr(primary.age, primary.brainage.data(:,1:12));
primary.mae = mean(abs(primary.age - primary.brainage.data(:,1:12)));
primary.wmae = mean(abs(primary.age - primary.brainage.data(:,1:12)))/primary.agerange;

secondary.age = secondary.demographics.data(:,2);
secondary.agerange = (max(secondary.age)-min(secondary.age));
secondary.rho = corr(secondary.age, secondary.brainage.data(:,1:12));
secondary.mae = mean(abs(secondary.age - secondary.brainage.data(:,1:12)));
secondary.wmae = mean(abs(secondary.age - secondary.brainage.data(:,1:12)))/secondary.agerange;

% set plot title
ttl{1,1} = 'gm xgtree';
ttl{1,2} = 'gm xglinear';
ttl{1,3} = 'gm rvm';
ttl{1,4} = 'gm stack';

ttl{2,1} = 'wm xgtree';
ttl{2,2} = 'wm xglinear';
ttl{2,3} = 'wm rvm';
ttl{2,4} = 'wm stack';

ttl{3,1} = 'gwm xgtree';
ttl{3,2} = 'gwm xglinear';
ttl{3,3} = 'gwm rvm';
ttl{3,4} = 'gwm stack';

% add correlations and mean absolute errors
primary.txt = cell(1);
secondary.txt = cell(1);
for i = 1:3
    for j = 1:4
        primary.txt{i,j} = compose(strcat({'MAE = '}, num2str(round(primary.mae(i*4-4+j),2),'%1.2f'), {' yrs\nwMAE = '}, num2str(round(primary.wmae(i*4-4+j),2),'%1.2f'), {' yrs\nrho = '},num2str(round(primary.rho(i*4-4+j),3),'%1.3f')));
        secondary.txt{i,j} = compose(strcat({'MAE = '}, num2str(round(secondary.mae(i*4-4+j),2),'%1.2f'), {' yrs\nwMAE = '}, num2str(round(secondary.wmae(i*4-4+j),2),'%1.2f'), {' yrs\nrho = '},num2str(round(secondary.rho(i*4-4+j),3),'%1.3f')));
    end
end
  
% make directory for figures
system('mkdir -p code/figures/');

% plot brain-predicted vs. chronological age (all models)
tissue = {'gm','wm','gwm'};

for i = 1:3
brainage_figure(i) = figure(); hold on;
for j = 1:4
    subplot(2,2,j); hold on;
    plot(primary.age, primary.brainage.data(:,i*4-4+j), '.','Color',[210/255, 210/255, 210/255],'Markersize',1,'LineWidth',0.25)
    plot(secondary.age, secondary.brainage.data(:,i*4-4+j), '.','Color',[0, 0.4470, 0.7410],'Markersize',4,'LineWidth',0.25)

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
    xticks(0:5:90);
    xticklabels(0:5:90);
    yticks(0:5:90);
    yticklabels(0:5:90);

    xh = get(gca,'xlabel');
    get(xh,'position');
    set(xh,'position', [63.5 39 -1]); 
    
    text(70,48,primary.txt{i,j},'Fontsize', 7, 'HorizontalAlignment', 'left');
    text(47,75,secondary.txt{i,j},'Fontsize', 7, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);
    title(ttl{i,j});
end

set(brainage_figure(i),'PaperUnits', 'centimeters');
set(brainage_figure(i),'PaperPosition', [0.25 7 20 20]);
print(strcat('code/figures/accuracy_', tissue{i}, '.pdf'),'-dpdf')
end
    
% plot brain-predicted vs. chronological age (stacked models)
ttl{1,4} = 'Grey matter';
ttl{2,4} = 'White matter';
ttl{3,4} = 'Grey and white matter';

j = 4;
brainage_figure(1) = figure(); hold on;
for i = 1:3
subplot(1,3,i); hold on;
    
    % scatter plot of brain-predicted vs. chronological age
    plot(primary.age, primary.brainage.data(:,i*4-4+j), '.','Color',[210/255, 210/255, 210/255],'Markersize',1,'LineWidth',0.25);
        xlim([min(primary.age) max(primary.age)]); ylim([min(primary.age) max(primary.age)]); 
        mdl = fitlm(primary.age,primary.brainage.data(:,i*4-4+j));
        prim_refline = refline(mdl.Coefficients{2,1}, mdl.Coefficients{1,1});
        prim_refline.Color = [180/255, 180/255, 180/255];
        
    plot(secondary.age, secondary.brainage.data(:,i*4-4+j), '.','Color',[0, 0.4470, 0.7410],'Markersize',2,'LineWidth',0.25);
        xlim([min(secondary.age) max(secondary.age)]); ylim([min(secondary.age) max(secondary.age)]); 
        mdl = fitlm(secondary.age,secondary.brainage.data(:,i*4-4+j));
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
    text(63.5,79.5,strcat(secondary.txt{i,j}, '{\hspace{0.4cm}}'), 'Fontsize', 7, 'HorizontalAlignment', 'right', 'Color', [0, 0.4470, 0.7410]);
    text(85,48,strcat(primary.txt{i,j}, '{\hspace{0.4cm}}'), 'Fontsize', 7, 'HorizontalAlignment', 'right');
    
    % add title
    title(ttl{i,j}, 'Fontsize', 9);
    
    % plot x=y identity line
    xy_identity = refline(1,0);
    xy_identity.Color = 'k';
end

set(brainage_figure(1),'PaperUnits', 'centimeters');
set(brainage_figure(1),'PaperSize',[20 6.5]);
set(brainage_figure(1),'PaperPosition', [0 0 20 6.5]);
print('code/figures/accuracy.pdf','-dpdf')
print('code/figures/accuracy2.png','-dpng', '-r300')
