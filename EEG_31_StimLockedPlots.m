% Fig S1
load([initPath.Exp '\data\group\eeg_stats\tfStat_trueUpVsShamUp.mat']) % from EEG_11_troughLockedAnalysis.m
load([initPath.Exp '\data\group\eeg_stats\tfStat_trueDownVsShamDown.mat']) % from EEG_11_troughLockedAnalysis.m
load([initPath.Exp 'data\group\eeg_stats\tfStat_trueUpvsTrueDown.mat']) % from EEG_11_troughLockedAnalysis.m
load([initPath.Exp '\data\group\eeg_stats\erpStat_trueUpVsShamUp.mat']) % from EEG_11_troughLockedAnalysis.m
load([initPath.Exp '\data\group\eeg_stats\erpStat_trueDownVsShamDown.mat']) % from EEG_11_troughLockedAnalysis.m
load([initPath.Exp 'data\group\eeg_stats\erpStat_trueUpvsTrueDown.mat']) % from EEG_11_troughLockedAnalysis.m

parulaHalfClear = [0.120903225806452,0.754387096774194,0.697670967741935;0.184429032258065,0.771745161290323,0.639109677419355;0.232335483870968,0.788816129032258,0.571925806451613;0.321229032258064,0.799632258064516,0.494625806451613;0.425522580645161,0.802896774193548,0.406574193548387;0.543361290322581,0.796035483870968,0.318687096774194;0.656267741935484,0.781867741935484,0.233203225806452;0.761322580645161,0.762416129032258,0.170558064516129;0.853903225806452,0.742577419354839,0.157364516129032;0.932683870967742,0.729283870967742,0.202977419354839;0.994006451612903,0.740158064516129,0.239851612903226;0.995622580645161,0.786170967741936,0.204903225806452;0.979764516129032,0.836200000000000,0.177732258064516;0.961296774193548,0.887400000000000,0.154329032258065;0.962651612903226,0.936567741935484,0.127070967741935;0.976900000000000,0.983900000000000,0.0805000000000000];
parulaHalfDark = [0.242200000000000,0.150400000000000,0.660300000000000;0.258006451612903,0.182303225806452,0.752561290322581;0.270987096774194,0.215877419354839,0.838800000000000;0.278574193548387,0.258025806451613,0.901500000000000;0.281390322580645,0.303467741935484,0.943674193548387;0.279009677419355,0.348196774193548,0.973687096774194;0.268354838709677,0.393622580645161,0.991625806451613;0.239206451612903,0.441225806451613,0.999390322580645;0.191209677419355,0.490777419354839,0.987087096774194;0.177783870967742,0.535067741935484,0.963974193548387;0.163893548387097,0.576790322580645,0.931422580645161;0.143777419354839,0.615677419354839,0.903645161290323;0.119487096774194,0.652822580645161,0.884212903225806;0.0858806451612903,0.686261290322581,0.851616129032258;0.0144677419354839,0.713783870967742,0.805416129032258;0.0197225806451613,0.736129032258065,0.753245161290323];
%%
cfg = [];
cfg.latency             = [-1.5 1.5];
cfg.frequency           = [5 30];
cfg.channel             = 'Cz';

tmpTftrueUp= ft_selectdata(cfg,grdAvg.trueUp);
tmpTftrueDown = ft_selectdata(cfg,grdAvg.trueDown);
tmpTfshamUp= ft_selectdata(cfg,grdAvg.shamUp);
tmpTfshamDown = ft_selectdata(cfg,grdAvg.shamDown);
tmpTfUpSub = ft_selectdata(cfg,grdAvg.subtractedUp);
tmpTfDownSub= ft_selectdata(cfg,grdAvg.subtractedDown);
tmpTfUpvsDown= ft_selectdata(cfg,grdAvg.subtractedUpvsDown);

tmpTfUpSub.mask = statTFUp.mask(2,:,:);
tmpTfDownSub.mask = statTFDown.mask(2,:,:);
tmpTfUpvsDown.mask = statTFupvsDown.mask(2,:,:);

cfg = [];
cfg.latency             = [-1.5 1.5];
cfg.channel             = 'Fz';
tmpErptrueUp= ft_selectdata(cfg,grdAvgERP.trueUp);
tmpErpshamUp= ft_selectdata(cfg,grdAvgERP.shamUp);
tmpErptrueDown= ft_selectdata(cfg,grdAvgERP.trueDown);
tmpErpshamDown= ft_selectdata(cfg,grdAvgERP.shamDown);
tmpErpsubtractedUp= ft_selectdata(cfg,grdAvgERP.substractedUp);
tmpErpsubtractedDown= ft_selectdata(cfg,grdAvgERP.substractedDown);



%% TF
cfg = [];
cfg.layout        = [initPath.FieldTrip '/template/layout/EEG1020.lay'];
cfg.maskstyle     = 'outline';
cfg.maskalpha     = 0.7;
cfg.fontsize      = 20;
cfg.zlim = [-0.1  0.4];
figure;ft_singleplotTFR(cfg,tmpTftrueUp); 
hold on ; yyaxis right; plot(tmpErptrueUp.time,tmpErptrueUp.avg,'-m','linewidth',2);
ylim([-60 40])
set(gca,'TickDir','both','TickLength',[0.015,0.015]); 
set(gcf,'color','w');
set(gca,'YColor',[0 0 0]);
box off 
xline(0,'-k','linewidth',1.5)
xline(0.1,'--k','linewidth',1.5)

figure;ft_singleplotTFR(cfg,tmpTftrueDown)
hold on ; yyaxis right; plot(tmpErptrueDown.time,tmpErptrueDown.avg,'-b','linewidth',2);
ylim([-60 40])
set(gca,'TickDir','both','TickLength',[0.015,0.015]); 
set(gcf,'color','w');
set(gca,'YColor',[0 0 0]);
box off 
xline(0,'-k','linewidth',1.5)
xline(0.1,'--k','linewidth',1.5)

figure;ft_singleplotTFR(cfg,tmpTfshamUp); 
hold on ; yyaxis right; plot(tmpErpshamUp.time,tmpErpshamUp.avg,'-y','linewidth',2);
ylim([-60 40])
set(gca,'TickDir','both','TickLength',[0.015,0.015]); 
set(gcf,'color','w');
set(gca,'YColor',[0 0 0]);
box off 
xline(0,'-k','linewidth',1.5)
xline(0.1,'--k','linewidth',1.5)

figure;ft_singleplotTFR(cfg,tmpTfshamDown)
hold on ; yyaxis right; plot(tmpErpshamDown.time,tmpErpshamDown.avg,'-c','linewidth',2);
ylim([-60 40])
set(gca,'TickDir','both','TickLength',[0.015,0.015]); 
set(gcf,'color','w');
set(gca,'YColor',[0 0 0]);
box off 
xline(0,'-k','linewidth',1.5)
xline(0.1,'--k','linewidth',1.5)


cfg.zlim = [-0.4  0.4];

cfg.maskparameter = 'mask';
figure;ft_singleplotTFR(cfg,tmpTfUpvsDown)
hold on ; yyaxis right; plot(tmpErptrueDown.time,tmpErptrueDown.avg,'-b','linewidth',2);
plot(tmpErptrueUp.time,tmpErptrueUp.avg,'-m','linewidth',2);
ylim([-60 40])
set(gca,'TickDir','both','TickLength',[0.015,0.015]); 
set(gcf,'color','w');
set(gca,'YColor',[0 0 0]);
box off 
xline(0,'-k','linewidth',1.5)
xline(0.1,'--k','linewidth',1.5)


cfg.zlim = [-0.2  0.2];

figure;ft_singleplotTFR(cfg,tmpTfUpSub); 
hold on ; yyaxis right; plot(tmpErptrueUp.time,tmpErptrueUp.avg,'-m','linewidth',2);
plot(tmpErpshamUp.time,tmpErpshamUp.avg,'-y','linewidth',2);
ylim([-60 40])
set(gca,'TickDir','both','TickLength',[0.015,0.015]); 
set(gcf,'color','w');
set(gca,'YColor',[0 0 0]);
box off 
xline(0,'-k','linewidth',1.5)
xline(0.1,'--k','linewidth',1.5)

figure;ft_singleplotTFR(cfg,tmpTfDownSub); 
hold on ; yyaxis right; plot(tmpErptrueDown.time,tmpErptrueDown.avg,'-b','linewidth',2);
plot(tmpErpshamDown.time,tmpErpshamDown.avg,'-c','linewidth',2);
ylim([-60 40])
set(gca,'TickDir','both','TickLength',[0.015,0.015]); 
set(gcf,'color','w');
set(gca,'YColor',[0 0 0]);
box off 
xline(0,'-k','linewidth',1.5)
xline(0.1,'--k','linewidth',1.5)

% 
% set(gca,'TickDir','both','TickLength',[0.015,0.015]); 
% set(gcf,'color','w');
% set(gca,'YColor',[0 0 0]);
% cbh = colorbar ; 
% cbh.Ticks = -1:0.5:1 ;
% cbh.TickDirection='both';
% cbh.TickLength=[0.015,0.015];
% box off 

%% SW Up 
cfg = [];
figure;ft_singleplotER(cfg,grdAvgERP.trueUp,grdAvgERP.trueDown,grdAvgERP.shamUp,grdAvgERP.shamDown)

figure;
%full SW
cfg = [];
cfg.channel= 'Fz';
cfg.xlim = [-1.5 1.5];
cfg.color = 'm';
h_plot_erf(cfg,allSubERP.trueUp');
cfg.color = 'y';
h_plot_erf(cfg,allSubERP.shamUp');
cfg.color = 'b';
h_plot_erf(cfg,allSubERP.trueDown');
cfg.color = 'c';
h_plot_erf(cfg,allSubERP.shamDown');

hold on
% Mark sign
tmpErptrueUp.nan= repmat(42,1,length(tmpErptrueUp.avg));
tmpErptrueUp.nan(find(staterpUp.mask(1,:)==0)) = NaN;
plot(tmpErptrueUp.time,tmpErptrueUp.nan,'-m','linewidth',2)

%Mark sign
tmpErptrueDown.nan= repmat(-63,1,length(tmpErptrueDown.avg));
tmpErptrueDown.nan(find(staterpDown.mask(1,:)==0)) = NaN;
plot(tmpErptrueDown.time,tmpErptrueDown.nan,'-b','linewidth',2)

%Mark sign
tmpErptrueDown.nan= repmat(-65,1,length(tmpErptrueDown.avg));
tmpErptrueDown.nan(find(staterpupvsDown.mask(1,:)==0)) = NaN;
plot(tmpErptrueDown.time,tmpErptrueDown.nan,'-k','linewidth',2)

%Plot characteristics
hold off
ylim([-67 45])
set(gca,'TickDir','both','TickLength',[0.015,0.015]); 
set(gcf,'color','w');
set(gca,'YColor',[0 0 0]);
box off 

xline(0,'-k','linewidth',1.5)
xline(0.1,'--k','linewidth',1.5)

