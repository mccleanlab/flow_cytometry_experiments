clearvars; clc; close all

%% Draw gate (comment out if loading previously saved data)
% gateParams.channels2gate = {'FSC-A', 'SSC-A'; 'FSC-A', 'FSC-H'};
% gateParams.channelsscale = {'linear','linear'; 'log','log'};
% gate = draw_gate(gateParams.channels2gate, gateParams.channelsscale);

%% Load previously saved gate
gate = load('gate_20210114_0_experiment_calibration_controls_yMM1454.mat');
gate = gate.gate_out;

%% Set plot options from file

% Get plot colors and Msn2_tex labels

% general_info_folder = 'D:\Google Drive\light_sweep_shared';
% general_info_folder = 'G:\My Drive\light_sweep_shared';
general_info_folder = 'D:\GoogleDriveUW\light_sweep_shared';
% general_info_folder = '/Volumes/GoogleDrive/My Drive/light_sweep_shared';
opts = detectImportOptions(fullfile(general_info_folder,'plot_settings.xlsx'),'Sheet','Msn2_tex');
opts.VariableTypes = {'categorical','categorical','double','double','double'};
plot_colors = readtable(fullfile(general_info_folder,'plot_settings.xlsx'),opts);
plot_colors.RGB = [plot_colors.R, plot_colors.G, plot_colors.B];
plot_colors.RGB(plot_colors.Msn2=='Msn2',:) = plot_colors.RGB(plot_colors.Msn2=='Msn2(WT|4E|WT)',:); % Replace Msn2 color
plot_colors = plot_colors(:,{'Msn2','Msn2_tex','RGB'});

%% Load previously calculated summary of data (comment out if running for 1st time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('data_stats.mat')
% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import measurements from .fcs files, label, and save to parquet files (uncomment if running for 1st time)
% folder_list = {
% %         fullfile('D:\Google Drive\flow_cytometry\20211129','plate_1')...
% %         fullfile('D:\Google Drive\flow_cytometry\20211129','plate_2')...
% %         fullfile('D:\Google Drive\flow_cytometry\20211202','plate_1')...
% %         fullfile('D:\Google Drive\flow_cytometry\20211202','plate_2')...
% %         fullfile('D:\Google Drive\flow_cytometry\20211207','plate_1')...
% %         fullfile('D:\Google Drive\flow_cytometry\20211207','plate_2')...
% 
%         fullfile('/Volumes/GoogleDrive/My Drive/flow_cytometry/20211129','plate_1')...
%         fullfile('/Volumes/GoogleDrive/My Drive/flow_cytometry/20211129','plate_2')...
%         fullfile('/Volumes/GoogleDrive/My Drive/flow_cytometry/20211202','plate_1')...
%         fullfile('/Volumes/GoogleDrive/My Drive/flow_cytometry/20211202','plate_2')...
%         fullfile('/Volumes/GoogleDrive/My Drive/flow_cytometry/20211207','plate_1')...
%         fullfile('/Volumes/GoogleDrive/My Drive/flow_cytometry/20211207','plate_2')...
%     };
% 
% for f = 1:numel(folder_list)
%     [data, experiment_name] = load_fcs('map','plate','folder',folder_list{f},'event_limit',25000,'plate_type','384_well');
%     data = add_gate(data,gate);
%     data = format_fcsdat(data);
%     parquetwrite(experiment_name,data);
% end
% 
% return

%% Load prevously saved data table from parquet files (comment out if running for 1st time)
data = import_parquet("*.parquet",{'experiment','well','FSCH','SSCH','BL2A','YL2A','strain','plasmid','Msn2','CLASP','reporter','condition','replicate','gate_net'});

%% Pre-process data

% Get rid of negative measurments
data(data.BL2A<=0,:) = [];
data(data.YL2A<=0,:) = [];

% Normalize mCitrine level by mScarlet level per cell
data.BL2A_norm = data.BL2A./data.YL2A;

% Change reporter naming to avoid tex issues
data.reporter = categorical(regexprep(string(data.reporter),'_',' '));

% Manually exclude outliers/weirdos
samples_to_exclude = [
    data.reporter=='glpT' & data.Msn2=='Msn2(S686A)' & data.replicate=='1',...
    data.reporter=='pALD3' & data.Msn2=='Msn2(S686A)' & data.replicate=='1',... % Abnormal mScarlet levels
    data.reporter=='pALD3' & data.Msn2=='Msn2(S686A)' & data.replicate=='2',... % Fails to induce in any condition
    data.reporter=='pTKL2' & data.Msn2=='Msn2(T,S686A)' & data.replicate=='4',...
    data.reporter=='pHXK1' & data.Msn2=='Msn2(T)' & data.replicate=='4',...
    ];
samples_to_exclude = any(samples_to_exclude,2);
data(samples_to_exclude,:) = []; % Note: only two replicates for ALD3/Msn2(S686A)

% Merge Msn2_tex labels into data
data = outerjoin(data, plot_colors(:,{'Msn2','Msn2_tex'}), 'Type', 'Left', 'MergeKeys', true);

%% Calculate population level statistics

% Calculate median expression and cell count per sample
data_stats = grpstats(data(data.gate_net==1,:),{'experiment','well','strain','reporter','plasmid','Msn2','Msn2_tex','CLASP','condition','replicate'},...
    {'nanmedian','nanmean','nanstd'},'DataVars',{'BL2A','YL2A'});
data_stats = clean_grpstats(data_stats,false);
data_stats.Properties.VariableNames('GroupCount') = {'cell_count'};
data_stats.Properties.VariableNames = regexprep(data_stats.Properties.VariableNames,'nanmedian_','');

% Merge cell count per well back into data
data = outerjoin(data, data_stats(:,{'experiment','strain','plasmid','condition','replicate','cell_count'}), 'Type', 'Left', 'MergeKeys', true);

% Calculate basal expression per sample
data_basal = data_stats(data_stats.condition=='control',{'experiment','strain','plasmid','replicate','BL2A','YL2A'});
data_basal.Properties.VariableNames(end-1:end) = {'BL2A_basal','YL2A_basal'};

% Automatically tag outliers based on basal mScarlet
[~,bound_low,bound_upper,~] = isoutlier(data_basal.YL2A_basal(data_basal.plasmid~='pMM0835'));
idx_keep = data_basal.YL2A_basal>bound_low & data_basal.YL2A_basal<bound_upper;
data_basal.exclude = ~idx_keep;
data_basal.exclude(data_basal.plasmid=='pMM0835') = 0; % Unexclude negative controls

% Merge indices for outliers back into main tables
data_stats = outerjoin(data_stats, data_basal, 'Type', 'Left', 'MergeKeys', true);
% data_stats.exclude(data_stats.cell_count<1000) = 1;
data = outerjoin(data, data_stats(:,{'experiment','strain','plasmid','replicate','BL2A_basal','YL2A_basal','exclude'}), 'Type', 'Left', 'MergeKeys', true);

% Calculate fold-change expression
data_stats.BL2A_fc = data_stats.BL2A./data_stats.BL2A_basal;
data_stats.YL2A_fc = data_stats.YL2A./data_stats.YL2A_basal;

% Calculate expression level for Msn2(S686A) per condition
grp_vars = {'strain','condition'};
data_Msn2_S686A = grpstats(data_stats(data_stats.Msn2=='Msn2(S686A)',:),grp_vars,'nanmean','DataVars','BL2A');
data_Msn2_S686A = clean_grpstats(data_Msn2_S686A);
data_Msn2_S686A.Properties.VariableNames(end) = {'BL2A_Msn2_S686A'};
data_stats = outerjoin(data_stats, data_Msn2_S686A, 'Type', 'Left', 'MergeKeys', true);

% Calculate expression noise
data_stats.BL2A_noise = (data_stats.nanstd_BL2A.^2)./(data_stats.nanmean_BL2A.^2);

condition_list = unique(data_stats.condition);

% Rename promoters to be less annoying
data_stats.reporter(data_stats.reporter~='glpT') = categorical(extractAfter(string(data_stats.reporter(data_stats.reporter~='glpT')),'p'));

%% Save data and summary table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('data.mat','data','-v7.3');
% save('data_stats.mat','data_stats')
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Figure 6: plot absolute Citrine all stresses (select reporters)
% reporters_to_plot = {'glpT','TKL2','ALD3','SIP18 D6','DCS2','HXK1','RTN2','SIP18','CTT1','SIP18 A4','DDR2','HSP12'};
reporters_to_plot = {'RTN2','DCS2','HSP12'};
Msn2_to_plot = {'Msn2(S686A)','Msn2','Msn2(T)','None'};
condition_set = {'control','0.5% glucose','0.1% glucose','0.01% glucose',...
    '0.0625 M NaCl','0.125 M NaCl','0.25 M NaCl','0.5 M NaCl'};
exclude = (data_stats.reporter=='ALD3' & data_stats.Msn2=='Msn2(S686A)'  & data_stats.replicate=='2') | ...
    (data_stats.reporter=='SIP18 A4' & data_stats.Msn2=='Msn2(S686A)'  & data_stats.replicate=='3');

clear g; close all; clc
figure('units','normalized','outerposition',[0 0 0.7 0.5]);
% figure('units','normalized','outerposition',[0 0 0.9 0.7]);
g = gramm('x',(data_stats.condition),'y',(data_stats.BL2A),'color',cellstr(data_stats.Msn2_tex),...
    'subset',~exclude & ismember(data_stats.Msn2,Msn2_to_plot) & data_stats.exclude==0 & data_stats.cell_count>1000);
g.facet_wrap(cellstr(data_stats.reporter),'ncols',6,'scale','independent');
g.stat_summary('type','std','geom',{'point','errorbar'},'setylim',true,'dodge',0.5,'width',2);
% g.set_names('x','','y','mCitrine','row','','column','','color','');
g.set_names('x','','y','','row','','column','','color','');
g.axe_property('XTickLabelRotation',45);
g.set_order_options('x',condition_set,'column',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_text_options('interpreter','tex');
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.geom_polygon('x',{[1.5 4.5]},'alpha',0.035);
g.no_legend();
g.draw();
g.redraw(0.05)
% g.facet_axes_handles(1).YLim = [0 750];
% export_fig(fullfile(pwd,'mCitrine_selected_reporters'),'-png','-m4');

%% Figure S12A: plot absolute Citrine all stresses (all reporters)
reporters_to_plot = {'glpT','TKL2','ALD3','SIP18 D6','DCS2','HXK1','RTN2','SIP18','CTT1','SIP18 A4','DDR2','HSP12'};
% reporters_to_plot = {'RTN2','TKL2','CTT1','HXK1','HSP12'};
% Msn2_to_plot = {'Msn2(S686A)','Msn2','Msn2(T)','None'};
Msn2_to_plot = unique(data_stats.Msn2);
condition_set = {'control','0.5% glucose','0.1% glucose','0.01% glucose',...
    '0.0625 M NaCl','0.125 M NaCl','0.25 M NaCl','0.5 M NaCl'};
exclude = (data_stats.reporter=='SIP18 A4' & data_stats.Msn2=='Msn2(S686A)'  & data_stats.replicate=='3');

clear g; close all; clc
figure('units','normalized','outerposition',[0 0 0.9 0.7]);
g = gramm('x',(data_stats.condition),'y',(data_stats.BL2A),'color',cellstr(data_stats.Msn2_tex),...
    'subset',~exclude & ismember(data_stats.Msn2,Msn2_to_plot) & data_stats.exclude==0 & data_stats.cell_count>1000);
g.facet_wrap(cellstr(data_stats.reporter),'ncols',6,'scale','independent');
g.stat_summary('type','std','geom',{'point','errorbar'},'setylim',true,'dodge',0.85,'width',3);
% g.set_names('x','','y','mCitrine','row','','column','','color','');
g.set_names('x','','y','','row','','column','','color','');
g.axe_property('XTickLabelRotation',45);
g.set_order_options('x',condition_set,'column',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_text_options('interpreter','tex');
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.geom_polygon('x',{[1.5 4.5]},'alpha',0.035);
% g.no_legend();
g.draw();
g.redraw(0.045)
g.facet_axes_handles(1).YLim = [0 750];
export_fig(fullfile(pwd,'mCitrine_all_reporters'),'-png','-m4');

%% Figure dissertation: plot absolute Citrine all stresses 
reporters_to_plot = {'glpT','TKL2','ALD3','SIP18 D6','DCS2','HXK1','RTN2','SIP18','CTT1','SIP18 A4','DDR2','HSP12'};
% reporters_to_plot = {'RTN2','DCS2','HSP12'};
Msn2_to_plot = {'Msn2(S686A)','Msn2','Msn2(S686D)','None'};
condition_set = {'control','0.01% glucose','0.5 M NaCl'};
exclude = (data_stats.reporter=='ALD3' & data_stats.Msn2=='Msn2(S686A)'  & data_stats.replicate=='2') | ...
    (data_stats.reporter=='SIP18 A4' & data_stats.Msn2=='Msn2(S686A)'  & data_stats.replicate=='3');
exclude = (data_stats.reporter=='SIP18 A4' & data_stats.Msn2=='Msn2(S686A)'  & data_stats.replicate=='3');

clear g; close all; clc
figure('units','normalized','outerposition',[0 0 0.7 0.5]);
% figure('units','normalized','outerposition',[0 0 0.9 0.7]);
g = gramm('x',(data_stats.condition),'y',(data_stats.BL2A),'color',cellstr(data_stats.Msn2_tex),...
    'subset',ismember(data_stats.Msn2,Msn2_to_plot) & ...
    ismember(data_stats.condition,condition_set) & data_stats.exclude==0 & data_stats.cell_count>1000);
g.facet_wrap(cellstr(data_stats.reporter),'ncols',6,'scale','independent');
g.stat_summary('type','std','geom',{'bar','black_errorbar'},'setylim',true,'dodge',0.75,'width',0.65);
g.set_names('x','','y','mCitrine','row','','column','','color','');
g.axe_property('XTickLabelRotation',45);
g.set_order_options('x',condition_set,'column',reporters_to_plot,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_text_options('interpreter','tex');
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.geom_polygon('x',{[1.5 4.5]},'alpha',0.035);
g.draw();
g.redraw(0.05)
g.facet_axes_handles(1).YLim = [0 750];
export_fig(fullfile(pwd,'absolute_mCitrine_S686_mutants'),'-png','-m4');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% OTHER PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot Msn2 levels per colony

clear g; close all; clc
figure('units','normalized','outerposition',[0 0 1 1]);
g = gramm('x',cellstr(data_stats.Msn2),'y',(data_stats.YL2A),...
    'color',double(data_stats.replicate),'label',(data_stats.replicate),...
    'marker',(data_stats.exclude | data_stats.cell_count<1000),...
    'subset',data_stats.condition=='control');
g.facet_wrap(cellstr(strcat(string(data_stats.strain),": ",string(data_stats.reporter))),'ncols',4,'scale','fixed');
g.geom_point('dodge',0.5);
g.set_names('x','','y','mScarlet','row','','column','','color','colony','marker','exclude');
g.axe_property('XTickLabelRotation',45,'TickLabelInterpreter','tex','YLim',[0 3500]);
g.set_order_options('x',0);
g.set_title('Msn2 levels per colony (20210707)');
g.set_point_options('markers',{'o','^'},'base_size',7)
g.draw();

g.update('y',(data_stats.YL2A) + 300);
g.geom_label('dodge',0.5);
g.no_legend();
g.draw();
g.export('file_name','Msn2_level_vs_colony','file_type','png');

%% Plot absolute mCitrine expression seperate the stresses 

% Plot absolute mCitrine (mean ± std)
Msn2_to_plot = {'Msn2(S686A)','Msn2','Msn2(T,S686A)','Msn2(T)','Msn2(S686D)','None',};
condition_sets  = {
    {'control','0.5% glucose','0.1% glucose','0.01% glucose'};...
    {'control','0.0625 M NaCl','0.125 M NaCl','0.25 M NaCl','0.5 M NaCl'}
    };
condition_labels = {'glucose','salt'};

for ii = 1:numel(condition_sets)
    condition = condition_sets{ii};
    label = condition_labels{ii};
    
    clear g; close all; clc
    figure('units','normalized','outerposition',[0 0 1 1]);
    g = gramm('x',(data_stats.Msn2),'y',log(data_stats.BL2A),'color',cellstr(data_stats.condition),...
        'subset',ismember(data_stats.Msn2,Msn2_to_plot) & ismember(data_stats.condition,condition)...
        & data_stats.exclude==0 & data_stats.cell_count>1000);
    g.facet_wrap(cellstr(data_stats.reporter),'ncols',6,'scale','independent');
    g.stat_summary('type','std','geom',{'point','errorbar'},'setylim',true);
    g.set_names('x','','y','log(mCitrine)','row','','column','','color','');
    g.axe_property('XTickLabelRotation',45);
    g.set_order_options('x',Msn2_to_plot,'color',condition);
    g.set_text_options('interpreter','tex');
    g.draw();
    
    for jj = 1:6
        g.facet_axes_handles(jj).XTickLabel = [];
    end
    
    g.redraw(0.025);
    %     export_fig(fullfile(pwd,strcat('meanstd_mCitrine_absolute_ZF_mutants_',label)),'-png','-m4');
end

%% Plot fold-change mCitrine seperate the stresses 

for ii = 1:numel(condition_sets)
    condition = condition_sets{ii};
    label = condition_labels{ii};
    
    clear g; close all; clc
    figure('units','normalized','outerposition',[0 0 1 1]);
    g = gramm('x',(data_stats.Msn2),'y',(data_stats.BL2A_fc),'color',cellstr(data_stats.condition),...
        'subset',ismember(data_stats.Msn2,Msn2_to_plot) & ismember(data_stats.condition,condition)...
        & data_stats.exclude==0 & data_stats.cell_count>1000);
    g.facet_wrap(cellstr(data_stats.reporter),'ncols',6,'scale','independent');
    g.stat_summary('type','std','geom',{'point','errorbar'},'setylim',true);
    g.set_names('x','','y','fold-change mCitrine','row','','column','','color','');
    g.axe_property('XTickLabelRotation',45);
    g.set_order_options('x',Msn2_to_plot,'color',condition);
    g.set_text_options('interpreter','tex');
    g.draw();
    
    for jj = 1:6
        g.facet_axes_handles(jj).XTickLabel = [];
    end
    
    g.redraw(0.025);
    export_fig(fullfile(pwd,strcat('meanstd_mCitrine_fc_ZF_mutants_',label)),'-png','-m4');
end

%% Plot fold change mCitrine all stresses
Msn2_to_plot = {'Msn2(S686A)','Msn2','Msn2(T,S686A)','Msn2(T)'};
condition_set = {'control','0.5% glucose','0.1% glucose','0.01% glucose',...
    '0.0625 M NaCl','0.125 M NaCl','0.25 M NaCl','0.5 M NaCl'};

clear g; close all; clc
figure('units','normalized','outerposition',[0 0 1 1]);
g = gramm('x',(data_stats.condition),'y',(data_stats.BL2A_fc),'color',cellstr(data_stats.Msn2_tex),...
    'subset',ismember(data_stats.Msn2,Msn2_to_plot) & data_stats.exclude==0 & data_stats.cell_count>1000);
g.facet_wrap(cellstr(data_stats.reporter),'ncols',4,'scale','independent');
g.stat_summary('type','std','geom',{'point','errorbar'},'setylim',true,'dodge',0.5,'width',1.5);
g.set_names('x','','y','fold-change mCitrine','row','','column','','color','');
g.axe_property('XTickLabelRotation',45);
g.set_order_options('x',condition_set,'color',plot_colors.Msn2_tex(ismember(plot_colors.Msn2,Msn2_to_plot)));
g.set_text_options('interpreter','tex');
g.set_color_options('map',plot_colors.RGB(ismember(plot_colors.Msn2,Msn2_to_plot),:)/255);
g.geom_polygon('x',{[1.5 4.5]},'alpha',0.035);
% g.geom_polygon('x',{[1.5 4.5]},'alpha',0.035);
g.draw();

for ii = 1:8
    g.facet_axes_handles(ii).XTickLabel = [];
end

for ii = 1:12
    g.facet_axes_handles(ii).YLim(1) = 0;
    if g.facet_axes_handles(ii).YLim(2) < 4
        g.facet_axes_handles(ii).YLim(2) = 4;
    end
end

g.redraw(0.025)

export_fig(fullfile(pwd,'mCitrine_fc_vs_stress_ZF_mutants'),'-png','-m4');



