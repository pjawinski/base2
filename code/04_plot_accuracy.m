% ==================================================
% === Plot brain-predicted vs. chronological age ===
% ==================================================
% /opt/matlab/bin/matlab -nodesktop -nodisplay 

% set working directory
cd /slow/projects/base2/

% load ukb and base2 data
primary = load('data/ukb_brainage.mat');
secondary = load('data/03_brainage.mat');

% calculate correlations and mean absolute errors
primary.age = primary.covs.data(:,5);
corr(primary.age, primary.brainage.data(:,1:12))
mean(abs(primary.age - primary.brainage.data(:,1:12)))

secondary.age = secondary.demographics.data(:,2);
corr(secondary.age, secondary.brainage.data(:,1:12))
mean(abs(secondary.age - secondary.brainage.data(:,1:12)))

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
txt = cell(1);
txt2 = cell(1);
for i = 1:3
    for j = 1:4
        txt{i,j} = compose(strcat({'MAE = '}, num2str(round(mean(abs(primary.age - primary.brainage.data(:,i*4-4+j))),2),'%1.2f'), {' yrs\n rho = '}, num2str(round(corr(primary.age, primary.brainage.data(:,i*4-4+j)),3),'%1.3f')));
        txt2{i,j} = compose(strcat({'MAE = '}, num2str(round(mean(abs(secondary.age - secondary.brainage.data(:,i*4-4+j))),2),'%1.2f'), {' yrs\n rho = '}, num2str(round(corr(secondary.age, secondary.brainage.data(:,i*4-4+j)),3),'%1.3f')));

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
    
    text(70,48,txt{i,j},'Fontsize', 7, 'HorizontalAlignment', 'left');
    text(47,75,txt2{i,j},'Fontsize', 7, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);
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

j = 4
brainage_figure(1) = figure(); hold on;
for i = 1:3
subplot(1,3,i); hold on;
    plot(primary.age, primary.brainage.data(:,i*4-4+j), '.','Color',[210/255, 210/255, 210/255],'Markersize',1,'LineWidth',0.25)
    plot(secondary.age, secondary.brainage.data(:,i*4-4+j), '.','Color',[0, 0.4470, 0.7410],'Markersize',2,'LineWidth',0.25)

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
    xticks(45:10:85);
    xticklabels(45:10:85);
    yticks(45:10:85);
    yticklabels(45:10:85);

    xh = get(gca,'xlabel');
    get(xh,'position');
    set(xh,'position', [63.5 37.5 -1]); 
    
    text(66,48,txt{i,j},'Fontsize', 7, 'HorizontalAlignment', 'left');
    text(45,78,txt2{i,j},'Fontsize', 7, 'HorizontalAlignment', 'left', 'Color', [0, 0.4470, 0.7410]);
    title(ttl{i,j}, 'Fontsize', 8);
    
end

set(brainage_figure(1),'PaperUnits', 'centimeters');
set(brainage_figure(1),'PaperPosition', [0.75 18 20 6]);
print('code/figures/accuracy.pdf','-dpdf')
