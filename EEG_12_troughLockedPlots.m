% Fig. 3 Fig. 4 Fig. S2

load([initPath.Exp '\data\group\eeg_stats\tfStat_UpVsDown.mat']) % from EEG_11_troughLockedAnalysis.m
load([initPath.Exp '\data\group\eeg_stats\tfStat_UpVsRest.mat']) % from EEG_11_troughLockedAnalysis.m
load([initPath.Exp '\data\group\eeg_stats\tfStat_DownVsRest.mat']) % from EEG_11_troughLockedAnalysis.m
load([initPath.Exp '\data\group\eeg_stats\erpStat_UpVsDown.mat']) % from EEG_11_troughLockedAnalysis.m
load([initPath.Exp '\data\group\eeg_stats\erpStat_UpVsRest.mat']) % from EEG_11_troughLockedAnalysis.m
load([initPath.Exp 'data\group\eeg_stats\erpStat_DownVsRest.mat']) % from EEG_11_troughLockedAnalysis.m
parulaHalfClear = [0.120903225806452,0.754387096774194,0.697670967741935;0.184429032258065,0.771745161290323,0.639109677419355;0.232335483870968,0.788816129032258,0.571925806451613;0.321229032258064,0.799632258064516,0.494625806451613;0.425522580645161,0.802896774193548,0.406574193548387;0.543361290322581,0.796035483870968,0.318687096774194;0.656267741935484,0.781867741935484,0.233203225806452;0.761322580645161,0.762416129032258,0.170558064516129;0.853903225806452,0.742577419354839,0.157364516129032;0.932683870967742,0.729283870967742,0.202977419354839;0.994006451612903,0.740158064516129,0.239851612903226;0.995622580645161,0.786170967741936,0.204903225806452;0.979764516129032,0.836200000000000,0.177732258064516;0.961296774193548,0.887400000000000,0.154329032258065;0.962651612903226,0.936567741935484,0.127070967741935;0.976900000000000,0.983900000000000,0.0805000000000000];
parulaHalfDark = [0.242200000000000,0.150400000000000,0.660300000000000;0.258006451612903,0.182303225806452,0.752561290322581;0.270987096774194,0.215877419354839,0.838800000000000;0.278574193548387,0.258025806451613,0.901500000000000;0.281390322580645,0.303467741935484,0.943674193548387;0.279009677419355,0.348196774193548,0.973687096774194;0.268354838709677,0.393622580645161,0.991625806451613;0.239206451612903,0.441225806451613,0.999390322580645;0.191209677419355,0.490777419354839,0.987087096774194;0.177783870967742,0.535067741935484,0.963974193548387;0.163893548387097,0.576790322580645,0.931422580645161;0.143777419354839,0.615677419354839,0.903645161290323;0.119487096774194,0.652822580645161,0.884212903225806;0.0858806451612903,0.686261290322581,0.851616129032258;0.0144677419354839,0.713783870967742,0.805416129032258;0.0197225806451613,0.736129032258065,0.753245161290323];

%%
cfg = [];
cfg.latency             = [-1.5 1.5];
cfg.frequency           = [5 30];
tmpTfUp= ft_selectdata(cfg,grdAvg_negPeak.up);
tmpTfDown = ft_selectdata(cfg,grdAvg_negPeak.down);
tmpTfRest = ft_selectdata(cfg,grdAvg_negPeak.rest);

cfg.channel             = 'Cz';

tmpTfUpVsDown = ft_selectdata(cfg,grdAvg_negPeak.upVsDown);
tmpTfUpVsRest = ft_selectdata(cfg,grdAvg_negPeak.upVsrest);
tmpTfDownVsRest = ft_selectdata(cfg,grdAvg_negPeak.downVsRest);

tmpTfUpVsDown.mask = statTFUpvsDown.mask(2,:,:);
tmpTfUpVsRest.mask = statTFUpvsRest.mask(2,:,:);
tmpTfDownVsRest.mask = statTFDownvsRest.mask(2,:,:);

cfg = [];
cfg.latency             = [-1.5 1.5];
cfg.channel             = 'Fz';
tmpErpUp= ft_selectdata(cfg,grdAvgERP_negPeak.up);
tmpErpDown= ft_selectdata(cfg,grdAvgERP_negPeak.down);
tmpErpRest = ft_selectdata(cfg,grdAvgERP_negPeak.rest);

%% TF/SW Up vs down
cfg = [];
cfg.layout        = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.maskstyle     = 'opacity';
cfg.maskalpha     = 0.6;
cfg.maskparameter = 'mask';
cfg.fontsize      = 20;
cfg.zlim = [-0.1 0.1];
figure;ft_singleplotTFR(cfg,tmpTfUpVsDown )


hold on
yyaxis right
%full SW
cfg = [];
cfg.channel= 'Fz';
cfg.xlim = [-1.5 1.5];
cfg.color = 'm';
h_plot_erf(cfg,allSubERP_negPeak.up');
cfg.color = 'b';
h_plot_erf(cfg,allSubERP_negPeak.down');

%Mark sign
tmpErpUp.nan= repmat(27,1,length(tmpErpUp.avg));
tmpErpUp.nan(find(staterpUpVsDown.mask(1,:)==0)) = NaN;
plot(tmpErpUp.time,tmpErpUp.nan,'-k','linewidth',2)


%Plot characteristics
hold off
ylim([-100 30])
set(gca,'TickDir','both','TickLength',[0.015,0.015]); 
set(gcf,'color','w');
set(gca,'YColor',[0 0 0]);
cbh = colorbar ; 
cbh.Ticks = -0.1:0.05:0.1 ;
cbh.TickDirection='both';
cbh.TickLength=[0.015,0.015];
box off 

% Up vs Rest
cfg = [];
cfg.layout        = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.maskstyle     = 'opacity';
cfg.maskalpha     = 0.6;
cfg.fontsize      = 20;
cfg.maskparameter = 'mask';
cfg.zlim = [-0.1 0.1];
figure;ft_singleplotTFR(cfg,tmpTfUpVsRest)


hold on
yyaxis right
%full SW
cfg = [];
cfg.channel= 'Fz';
cfg.xlim = [-1.5 1.5];
cfg.color = 'm';
h_plot_erf(cfg,allSubERP_negPeak.up');
cfg.color = 'g';
h_plot_erf(cfg,allSubERP_negPeak.rest');


%Mark sign
tmpErpRest.nan= repmat(27,1,length(tmpErpRest.avg));
tmpErpRest.nan(find(staterpUpVsRest.mask(1,:)==0)) = NaN;
plot(tmpErpRest.time,tmpErpRest.nan,'-k','linewidth',2)



%Plot characteristics
hold off
ylim([-100 30])
set(gca,'TickDir','both','TickLength',[0.015,0.015]); 
set(gcf,'color','w');
set(gca,'YColor',[0 0 0]);
cbh = colorbar ; 
cbh.Ticks = -0.1:0.05:0.1 ;
cbh.TickDirection='both';
cbh.TickLength=[0.015,0.015];
box off 

% Down vs Rest
cfg = [];
cfg.layout        = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.maskstyle     = 'opacity';
cfg.maskalpha     = 0.6;
cfg.maskparameter = 'mask';
cfg.fontsize      = 20;
cfg.zlim = [-0.1 0.1];
figure;ft_singleplotTFR(cfg,tmpTfDownVsRest)


hold on
yyaxis right
%full SW
cfg = [];
cfg.channel= 'Fz';
cfg.xlim = [-1.5 1.5];
cfg.color = 'b';
h_plot_erf(cfg,allSubERP_negPeak.down');
cfg.color = 'g';
h_plot_erf(cfg,allSubERP_negPeak.rest');


%Mark sign
tmpErpRest.nan= repmat(27,1,length(tmpErpRest.avg));
tmpErpRest.nan(find(staterpDownVsRest.mask(1,:)==0)) = NaN;
plot(tmpErpRest.time,tmpErpRest.nan,'-k','linewidth',2)

%Plot characteristics
hold off
ylim([-100 30])
set(gca,'TickDir','both','TickLength',[0.015,0.015]); 
set(gcf,'color','w');
set(gca,'YColor',[0 0 0]);
cbh = colorbar ; 
cbh.Ticks = -0.1:0.05:0.1 ;
cbh.TickDirection='both';
cbh.TickLength=[0.015,0.015];
box off 


%% Topoplot TF Up vs Down
cfg = [];
cfg.xlim = [0.25 0.4];
cfg.ylim = [12 17];
cfg.zlim = [0 0.08];
cfg.colormap = parulaHalfClear;
cfg.layout = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.colorbar           = 'yes';
cfg.style              = 'straight';
figure; cfg.title = 'up vs down';ft_topoplotTFR(cfg,grdAvg_negPeak.upVsDown); 
set(gcf,'color','w');

% Topoplot TF Up vs rest
cfg = [];
cfg.xlim = [-0.25 0.08];
cfg.ylim = [7 12];
cfg.zlim = [-0.15 0];
cfg.colormap = parulaHalfDark;
cfg.layout = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.colorbar           = 'yes';
cfg.style              = 'straight';
figure; cfg.title = 'up vs rest'; ft_topoplotTFR(cfg,grdAvg_negPeak.upVsrest); 
set(gcf,'color','w');

% Topoplot TF down vs rest (neg cluster)
cfg = [];
cfg.xlim = [-0.2 0.01];
cfg.ylim = [8 1];
cfg.zlim = [-0.15 0];
cfg.colormap = parulaHalfDark;
cfg.layout = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.colorbar           = 'yes';
cfg.style              = 'straight';
figure; cfg.title = 'down vs rest neg'; ft_topoplotTFR(cfg,grdAvg_negPeak.downVsRest); 
set(gcf,'color','w');

% Topoplot TF down vs rest (pos cluster)
cfg = [];
cfg.xlim = [0.8 1.25];
cfg.ylim = [5 8];
cfg.zlim = [0 0.05];
cfg.colormap = parulaHalfClear;
cfg.layout = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.colorbar           = 'yes';
cfg.style              = 'straight';
figure; cfg.title = 'down vs rest pos'; ft_topoplotTFR(cfg,grdAvg_negPeak.downVsRest); 
set(gcf,'color','w');


%% Zoom SW 

%Up vs Down
figure ; cfg = [];
cfg.channel= 'Fz';
cfg.xlim = [0.2 1.5];
cfg.color = 'm';
h_plot_erf(cfg,allSubERP_negPeak.up');
cfg.color = 'b';
h_plot_erf(cfg,allSubERP_negPeak.down');
hold on 
% %Mark sign
tmpErpUp.nan= tmpErpUp.avg;
tmpErpDown.nan= tmpErpDown.avg;
tmpErpUp.nan(find(staterpUpVsDown.mask(1,:)==0)) = NaN;
tmpErpDown.nan(find(staterpUpVsDown.mask(1,:)==0)) = NaN;
plot(tmpErpUp.time,tmpErpUp.nan,'-m','linewidth',3)
plot(tmpErpDown.time,tmpErpDown.nan,'-b','linewidth',3)

%Plot characteristics
hold off
ylim([-15 25])
set(gca,'TickDir','both','TickLength',[0.015,0.015]); 
set(gcf,'color','w');
box off 


% Up vs rest
figure ; cfg = [];
cfg.channel= 'Fz';%allSubERP.all{1}.label{idx_channel};
cfg.xlim = [0.2 1.5];
cfg.color = 'm';
h_plot_erf(cfg,allSubERP_negPeak.up');
cfg.color = 'g';
h_plot_erf(cfg,allSubERP_negPeak.rest');
hold on 
% %Mark sign
tmpErpUp.nan= tmpErpUp.avg;
tmpErpRest.nan= tmpErpRest.avg;
tmpErpUp.nan(find(staterpUpVsRest.mask(1,:)==0)) = NaN;
tmpErpRest.nan(find(staterpUpVsRest.mask(1,:)==0)) = NaN;
plot(tmpErpUp.time,tmpErpUp.nan,'-m','linewidth',3)
plot(tmpErpRest.time,tmpErpRest.nan,'-g','linewidth',3)

%Plot characteristics
hold off
ylim([-15 25])
set(gca,'TickDir','both','TickLength',[0.015,0.015]); 
set(gcf,'color','w');
box off 

% Down  vs rest
figure ; cfg = [];
cfg.channel= 'Fz';
cfg.xlim = [0.2 1.5];
cfg.color = 'b';
h_plot_erf(cfg,allSubERP_negPeak.down');
cfg.color = 'g';
h_plot_erf(cfg,allSubERP_negPeak.rest');
hold on 
% %Mark sign
tmpErpDown.nan= tmpErpDown.avg;
tmpErpRest.nan= tmpErpRest.avg;
tmpErpDown.nan(find(staterpDownVsRest.mask(1,:)==0)) = NaN;
tmpErpRest.nan(find(staterpDownVsRest.mask(1,:)==0)) = NaN;
plot(tmpErpDown.time,tmpErpDown.nan,'-b','linewidth',3)
plot(tmpErpRest.time,tmpErpRest.nan,'-g','linewidth',3)

%Plot characteristics
hold off
ylim([-15 25])
set(gca,'TickDir','both','TickLength',[0.015,0.015]); 
set(gcf,'color','w');
box off 

%% Topoplot ERP  Peak

%Up vs Down and Up vs Rest and Down vs Rest

cfg = [];
cfg.xlim     =  [0.32 0.64];
cfg.zlim     = [0 20];
cfg.style    = 'straight';
cfg.colorbar = 'yes';
cfg.colormap = parulaHalfClear;
cfg.layout   = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
figure; cfg.title = 'up'; ft_topoplotER(cfg,grdAvgERP_negPeak.up); set(gcf,'color','w');
figure; cfg.title = 'down'; ft_topoplotER(cfg,grdAvgERP_negPeak.down); set(gcf,'color','w');
figure; cfg.title = 'rest'; ft_topoplotER(cfg,grdAvgERP_negPeak.rest); set(gcf,'color','w');

cfg.xlim     =  [0.32 0.64];
cfg.zlim     = [0 3];
figure; cfg.title = 'up vs down';ft_topoplotER(cfg,grdAvgERP_negPeak.upVsDown); set(gcf,'color','w');
figure; cfg.title = 'up vs rest';ft_topoplotER(cfg,grdAvgERP_negPeak.upVsrest); set(gcf,'color','w');
figure; cfg.title = 'down vs rest';ft_topoplotER(cfg,grdAvgERP_negPeak.downVsRest); set(gcf,'color','w');

%% Topoplot ERP  second cluster

% Up vs Down 

cfg = [];
cfg.xlim     =  [0.78 1.03];
cfg.zlim     = [-2 2];
cfg.style    = 'straight';
cfg.colorbar = 'yes';
cfg.layout   = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.colormap = parula;
figure; cfg.title = 'up'; ft_topoplotER(cfg,grdAvgERP_negPeak.up); set(gcf,'color','w');
figure; cfg.title = 'down'; ft_topoplotER(cfg,grdAvgERP_negPeak.down); set(gcf,'color','w');

cfg.zlim     = [-3 0];
cfg.colormap = parulaHalfDark;
figure; cfg.title = 'up vs down'; ft_topoplotER(cfg,grdAvgERP_negPeak.upVsDown); set(gcf,'color','w');

% Up vs Rest
cfg = [];
cfg.xlim     =  [0.72 1.14];
cfg.zlim     = [-2 2];
cfg.style    = 'straight';
cfg.colorbar = 'yes';
cfg.layout   = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.colormap = parula;
figure; cfg.title = 'up'; ft_topoplotER(cfg,grdAvgERP_negPeak.up); set(gcf,'color','w');
figure; cfg.title = 'rest'; ft_topoplotER(cfg,grdAvgERP_negPeak.rest); set(gcf,'color','w');

cfg.zlim     = [-3 0];
cfg.colormap = parulaHalfDark;
figure; cfg.title = 'up vs rest';ft_topoplotER(cfg,grdAvgERP_negPeak.upVsrest); set(gcf,'color','w');

% Down vs Rest
cfg = [];
cfg.xlim     =  [1.18 1.5];
cfg.zlim     = [0 3];
cfg.style    = 'straight';
cfg.colorbar = 'yes';
cfg.layout   = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.colormap = parulaHalfClear;
figure; cfg.title = 'down'; ft_topoplotER(cfg,grdAvgERP_negPeak.down); set(gcf,'color','w');
figure; cfg.title = 'rest'; ft_topoplotER(cfg,grdAvgERP_negPeak.rest); set(gcf,'color','w');

cfg.zlim     = [-3 0];
cfg.colormap = parulaHalfDark;
figure; cfg.title = 'down vs rest';ft_topoplotER(cfg,grdAvgERP_negPeak.downVsRest); set(gcf,'color','w');




%% TF/SW Up vs down

for idx_cond = 1 : 3
    for idx_chan = 1: 6
        
        cfg = [];
        cfg.layout        = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
        cfg.fontsize      = 20;
        cfg.zlim = [-0.5 0.5];
        cfg.channel= tmpTfUp.label{idx_chan};

        if idx_cond == 1
            cfg.title = ['up ' tmpTfUp.label{idx_chan}];
            figure;ft_singleplotTFR(cfg,tmpTfUp)
        elseif idx_cond == 2
            cfg.title = ['down ' tmpTfUp.label{idx_chan}];
            figure;ft_singleplotTFR(cfg,tmpTfDown)
        elseif idx_cond == 3
            cfg.title = ['rest ' tmpTfUp.label{idx_chan}];
            figure;ft_singleplotTFR(cfg,tmpTfRest)
        end

        hold on
        yyaxis right
        %full SW
        cfg = [];
        cfg.channel= tmpTfUp.label{idx_chan};
        cfg.xlim = [-1.5 1.5];
        if idx_cond == 1
            cfg.color = 'm';
            h_plot_erf(cfg,allSubERP_negPeak.up');
        elseif idx_cond == 2
            cfg.color = 'b';
            h_plot_erf(cfg,allSubERP_negPeak.down');
        elseif idx_cond == 3
            cfg.color = 'g';
            h_plot_erf(cfg,allSubERP_negPeak.rest');
        end
                                
        %Plot characteristics
        hold off
        ylim([-100 30])
        set(gca,'TickDir','both','TickLength',[0.015,0.015]);
        set(gcf,'color','w');
        set(gca,'YColor',[0 0 0]);
        cbh = colorbar ;
%         cbh.Ticks = -0.1:0.05:0.1 ;
%         cbh.TickDirection='both';
%         cbh.TickLength=[0.015,0.015];
        box off
        
    end
end


%% Per channel stats


cfg = [];
cfg.latency             = [-1.5 1.5];
cfg.frequency           = [5 30];
cfr.avgoverchan = 'no';
tmpTfUpVsDown = ft_selectdata(cfg,grdAvg_negPeak.upVsDown);
tmpTfUpVsRest = ft_selectdata(cfg,grdAvg_negPeak.upVsrest);
tmpTfDownVsRest = ft_selectdata(cfg,grdAvg_negPeak.downVsRest);

tmpTfUpVsDown.mask = statTFUpvsDown.mask;
tmpTfDownVsRest.mask = statTFDownvsRest.mask;
tmpTfUpVsRest.mask = statTFUpvsRest.mask;


for idx_cond = 3%1 : 3
    for idx_chan = 1%1   : 6
        
        cfg = [];
        cfg.layout        = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
        cfg.fontsize      = 20;
        cfg.zlim = [-0.5 0.5];
        cfg.channel= tmpTfUpVsDown.label{idx_chan};
        cfg.maskalpha     = 0.6;
        cfg.maskparameter = 'mask';

        if idx_cond == 1
            cfg.title = ['up vs down' tmpTfUpVsDown.label{idx_chan}];
            figure;ft_singleplotTFR(cfg,tmpTfUpVsDown)
        elseif idx_cond == 2
            cfg.title = ['up vs rest ' tmpTfUpVsDown.label{idx_chan}];
            figure;ft_singleplotTFR(cfg,tmpTfUpVsRest)
        elseif idx_cond == 3
            cfg.title = ['down vs rest ' tmpTfUpVsDown.label{idx_chan}];
            figure;ft_singleplotTFR(cfg,tmpTfDownVsRest)
        end

        

        hold on
        yyaxis right
        
        %full SW
        cfg = [];
        cfg.channel= tmpTfUpVsDown.label{idx_chan};
        cfg.xlim = [-1.5 1.5];
        if idx_cond == 1
            cfg.color = 'm';
            h_plot_erf(cfg,allSubERP_negPeak.up');
            cfg.color = 'b';
            h_plot_erf(cfg,allSubERP_negPeak.down');
            %Mark sign
            tmpErpRest.nan= repmat(27,1,length(tmpErpRest.avg));
            tmpErpRest.nan(find(staterpUpVsDown.mask(idx_chan,:)==0)) = NaN;
            plot(tmpErpRest.time,tmpErpRest.nan,'-k','linewidth',2)
        elseif idx_cond == 2
            cfg.color = 'm';
            h_plot_erf(cfg,allSubERP_negPeak.up');
            cfg.color = 'g';
            h_plot_erf(cfg,allSubERP_negPeak.rest');
            tmpErpRest.nan= repmat(27,1,length(tmpErpRest.avg));
            tmpErpRest.nan(find(staterpUpVsRest.mask(idx_chan,:)==0)) = NaN;
            plot(tmpErpRest.time,tmpErpRest.nan,'-k','linewidth',2)

        elseif idx_cond == 3
            cfg.color = 'b';
            h_plot_erf(cfg,allSubERP_negPeak.down');
            cfg.color = 'g';
            h_plot_erf(cfg,allSubERP_negPeak.rest');
            tmpErpRest.nan= repmat(27,1,length(tmpErpRest.avg));
            tmpErpRest.nan(find(staterpDownVsRest.mask(idx_chan,:)==0)) = NaN;
            plot(tmpErpRest.time,tmpErpRest.nan,'-k','linewidth',2)

        end
                                
        %Plot characteristics
        hold off
        ylim([-100 30])
        set(gca,'TickDir','both','TickLength',[0.015,0.015]);
        set(gcf,'color','w');
        set(gca,'YColor',[0 0 0]);
        cbh = colorbar ;
        box off
        hold off
    end
end


