%% Load data & baseline
clc
[listSub,listSubEEG,listSubBehav] = getFullDatasets;
%Utilities
parulaHalfClear = [0.120903225806452,0.754387096774194,0.697670967741935;0.184429032258065,0.771745161290323,0.639109677419355;0.232335483870968,0.788816129032258,0.571925806451613;0.321229032258064,0.799632258064516,0.494625806451613;0.425522580645161,0.802896774193548,0.406574193548387;0.543361290322581,0.796035483870968,0.318687096774194;0.656267741935484,0.781867741935484,0.233203225806452;0.761322580645161,0.762416129032258,0.170558064516129;0.853903225806452,0.742577419354839,0.157364516129032;0.932683870967742,0.729283870967742,0.202977419354839;0.994006451612903,0.740158064516129,0.239851612903226;0.995622580645161,0.786170967741936,0.204903225806452;0.979764516129032,0.836200000000000,0.177732258064516;0.961296774193548,0.887400000000000,0.154329032258065;0.962651612903226,0.936567741935484,0.127070967741935;0.976900000000000,0.983900000000000,0.0805000000000000];


%%
allSubERP_negPeak = {};
allSub_negPeak = {};
                                                           
for idx_sub = 1 :  length(listSubEEG)
    
    sub = listSubEEG{idx_sub};
    
    disp (['loading ' sub ' dataset'])
    
    load ([initPath.Exp '\data\' sub '\exp\eeg\' sub '_preprocessed_continuous.mat'])
    
    
    trl = load ([initPath.Exp '\data\' sub '\exp\eeg\' sub '_trl_trough_EEG.csv']); % time of negPeak / 1 = down ; 2 = rest; 3 = up; / sleep stage
    
    cfg=[];
    cfg.trl    = [trl(:,1)-2*data.fsample trl(:,1)+2*data.fsample repmat(-2*data.fsample,length(trl(:,1)),1)];
    data_epoched = ft_redefinetrial(cfg, data);
    
    cfg = [];
    cfg.resamplefs      = 100;
    cfg.demean          = 'yes';
    data_epoched = ft_resampledata(cfg, data_epoched);
    
    data_epoched.trialinfo = trl(:,2:3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % ERP
    % NREM up vs down and all
    
    %down
    cfg = [];
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,1)== 1 & (data_epoched.trialinfo(:,2)== 2 | data_epoched.trialinfo(:,2)== 3));
    timelock_data= ft_timelockanalysis(cfg, data_epoched);
    timelock_data= rmfield(timelock_data,'cfg');
    
    allSubERP_negPeak.nbTrl(idx_sub,2) =  length(cfg.trials);
    timelock_data.cfg=[];
    allSubERP_negPeak.down{idx_sub}= timelock_data;
    
    %rest
    
    cfg = [];
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,1)== 2 & (data_epoched.trialinfo(:,2)== 2 | data_epoched.trialinfo(:,2)== 3));
    timelock_data= ft_timelockanalysis(cfg, data_epoched);
    timelock_data= rmfield(timelock_data,'cfg');
    
    allSubERP_negPeak.nbTrl(idx_sub,3) =  length(cfg.trials);
    timelock_data.cfg=[];
    allSubERP_negPeak.rest{idx_sub}= timelock_data;
    
    %up
    cfg = [];
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,1)== 3 & (data_epoched.trialinfo(:,2)== 2 | data_epoched.trialinfo(:,2)== 3));
    timelock_data= ft_timelockanalysis(cfg, data_epoched);
    timelock_data= rmfield(timelock_data,'cfg');
    
    allSubERP_negPeak.nbTrl(idx_sub,1) =  length(cfg.trials);
    timelock_data.cfg=[];
    allSubERP_negPeak.up{idx_sub}= timelock_data;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % TF
    
    % NREM up vs down and all

    foi = 5:0.5:30; % 0.1 Hz steps
    toi = -3:0.01:3; % 0.1 s steps

    %down
    cfg = [];
    cfg.method     = 'mtmconvol';
    cfg.pad        = 'nextpow2';
    cfg.taper      = 'hanning';
    cfg.foi        = foi; % 0.1 Hz steps
    cfg.toi        = toi; % 0.1 s steps
    cfg.t_ftimwin  = 5./cfg.foi;
    cfg.tapsmofrq  = 0.4 *cfg.foi;
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,1)== 1 & (data_epoched.trialinfo(:,2)== 2 | data_epoched.trialinfo(:,2)== 3));
    timelock_data  = ft_freqanalysis(cfg, data_epoched);

    allSub_negPeak.nbTrl(idx_sub,2) =  length(cfg.trials);

    cfg = [];
    cfg.baseline       = [-2 2];
    cfg.baselinetype = 'relchange';
    timelock_data    = ft_freqbaseline(cfg, timelock_data);
    timelock_data= rmfield(timelock_data,'cfg');

    timelock_data.cfg        =[];
    timelock_data.freq       = foi;
    allSub_negPeak.down{idx_sub}= timelock_data;

    %rest
    cfg = [];
    cfg.method     = 'mtmconvol';
    cfg.pad        = 'nextpow2';
    cfg.taper      = 'hanning';
    cfg.foi        = foi; % 0.1 Hz steps
    cfg.toi        = toi; % 0.1 s steps
    cfg.t_ftimwin  = 5./cfg.foi;
    cfg.tapsmofrq  = 0.4 *cfg.foi;
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,1)== 2 & (data_epoched.trialinfo(:,2)== 2 | data_epoched.trialinfo(:,2)== 3));
    timelock_data  = ft_freqanalysis(cfg, data_epoched);

    allSub_negPeak.nbTrl(idx_sub,3) =  length(cfg.trials);

    cfg = [];
    cfg.baseline       = [-2 2];
    cfg.baselinetype = 'relchange';
    timelock_data    = ft_freqbaseline(cfg, timelock_data);
    timelock_data= rmfield(timelock_data,'cfg');

    timelock_data.cfg        =[];
    timelock_data.freq       = foi;
    allSub_negPeak.rest{idx_sub}= timelock_data;

        %up
    cfg = [];
    cfg.method     = 'mtmconvol';
    cfg.pad        = 'nextpow2';
    cfg.taper      = 'hanning';
    cfg.foi        = foi;
    cfg.toi        = toi;
    cfg.t_ftimwin  = 5./cfg.foi;
    cfg.tapsmofrq  = 0.4 *cfg.foi;
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,1)== 3 & (data_epoched.trialinfo(:,2)== 2 | data_epoched.trialinfo(:,2)== 3));
    timelock_data  = ft_freqanalysis(cfg, data_epoched);

    allSub_negPeak.nbTrl(idx_sub,1) =  length(cfg.trials);

    cfg = [];
    cfg.baseline       = [-2 2];
    cfg.baselinetype = 'relchange';
    timelock_data    = ft_freqbaseline(cfg, timelock_data);
    timelock_data= rmfield(timelock_data,'cfg');

    timelock_data.cfg        =[];
    timelock_data.freq       = foi;
    allSub_negPeak.up{idx_sub}= timelock_data;


    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'x1-x2';
    allSub_negPeak.upVsDown{idx_sub}   = ft_math(cfg, allSub_negPeak.up{idx_sub},allSub_negPeak.down{idx_sub});
    allSub_negPeak.upVsrest{idx_sub}   = ft_math(cfg, allSub_negPeak.up{idx_sub},allSub_negPeak.rest{idx_sub});
    allSub_negPeak.downVsRest{idx_sub} = ft_math(cfg, allSub_negPeak.down{idx_sub},allSub_negPeak.rest{idx_sub});
end



cfg = [];
cfg.keepindividual = 'yes'; % for cohen's d computation
grdAvgERP_negPeak.up = ft_timelockgrandaverage(cfg, allSubERP_negPeak.up{:});
grdAvgERP_negPeak.down = ft_timelockgrandaverage(cfg, allSubERP_negPeak.down{:});
grdAvgERP_negPeak.rest = ft_timelockgrandaverage(cfg, allSubERP_negPeak.rest{:});
cfg = [];
cfg.operation='x1-x2';
cfg.parameter='avg';
grdAvgERP_negPeak.upVsDown = ft_math(cfg, grdAvgERP_negPeak.up, grdAvgERP_negPeak.down);
grdAvgERP_negPeak.upVsrest = ft_math(cfg, grdAvgERP_negPeak.up, grdAvgERP_negPeak.rest);
grdAvgERP_negPeak.downVsRest = ft_math(cfg, grdAvgERP_negPeak.down, grdAvgERP_negPeak.rest);


cfg = [];
cfg.keepindividual = 'yes'; % for cohen's d computation
grdAvg_negPeak.up = ft_freqgrandaverage(cfg, allSub_negPeak.up{:});
grdAvg_negPeak.down = ft_freqgrandaverage(cfg, allSub_negPeak.down{:});
grdAvg_negPeak.rest = ft_freqgrandaverage(cfg, allSub_negPeak.rest{:});
grdAvg_negPeak.upVsDown = ft_freqgrandaverage(cfg, allSub_negPeak.upVsDown{:});
grdAvg_negPeak.upVsrest = ft_freqgrandaverage(cfg, allSub_negPeak.upVsrest{:});
grdAvg_negPeak.downVsRest = ft_freqgrandaverage(cfg, allSub_negPeak.downVsRest{:});

csvwrite([initPath.Exp '\data\group\nbTrl_anlayzed.CSV'],allSubERP_negPeak.nbTrl)



%% stat ERP

load neighboursPerso.mat

cfg                     = [];
cfg.design(1,1:2*(length(listSubEEG)))  = [ones(1,(length(listSubEEG))) 2*ones(1,(length(listSubEEG)))];
cfg.design(2,1:2*(length(listSubEEG)))  = [1:(length(listSubEEG)) 1:(length(listSubEEG))];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
cfg.method              = 'montecarlo';       
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'cluster';
cfg.clusterstatistic    = 'maxsum'; 
cfg.minnbchan           = 0;              
cfg.neighbours          = neighbours_perso; 
cfg.tail                = 0;                    
cfg.clustertail         = 0;
cfg.alpha               = 0.025/3;%0.05/6; 
cfg.clusteralpha        = 0.01; 
% cfg.tfce_H              = 2;       % default setting
% cfg.tfce_E              = 0.5;     % default setting
cfg.numrandomization    = 500;     
cfg.latency             = [-1.5 1.5];


[staterpUpVsDown] = ft_timelockstatistics(cfg,  allSubERP_negPeak.up{:}, allSubERP_negPeak.down{:});
[staterpUpVsRest] = ft_timelockstatistics(cfg,  allSubERP_negPeak.up{:}, allSubERP_negPeak.rest{:});
[staterpDownVsRest] = ft_timelockstatistics(cfg,  allSubERP_negPeak.down{:}, allSubERP_negPeak.rest{:});

save([initPath.Exp '\data\group\eeg_stats\erpStat_UpVsDown.mat'],staterpUpVsDown) % for EEG12
save([initPath.Exp '\data\group\eeg_stats\erpStat_UpVsRest.mat'],staterpUpVsRest) % for EEG12
save([initPath.Exp '\data\group\eeg_stats\erpStat_DownVsRest.mat'],staterpDownVsRest) % for EEG12

% Peak amplitude  for correlation (IRM)
peakAmplitude = [];

for idx_sub = 1 : length(allSubERP_negPeak.up)

    cfg = [];
    cfg.avgoverchan         = 'yes';
    cfg.latency             = [0.32 0.64];
    cfg.avgovertime         = 'yes';

    tmp = ft_selectdata(cfg,allSubERP_negPeak.up{idx_sub});

    peakAmplitude(idx_sub,1) = mean(tmp.avg);
    
    tmp = ft_selectdata(cfg,allSubERP_negPeak.down{idx_sub});

    peakAmplitude(idx_sub,2) = mean(tmp.avg);

  
end

csvwrite([initPath.Exp '\data\group\peakAmplitudeContrastUpvsDown.csv'],peakAmplitude) %=> to write getEEG() function

%% Stat TF

load neighboursPerso.mat

cfg                     = [];
cfg.design(1,1:2*length(listSubEEG))  = [ones(1,length(listSubEEG)) 2*ones(1,length(listSubEEG))];
cfg.design(2,1:2*length(listSubEEG))  = [1:length(listSubEEG) 1:length(listSubEEG)];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'cluster';
cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 0;
cfg.neighbours          = neighbours_perso;
cfg.alpha               = 0.025/3;
cfg.clusteralpha        = 0.01;
cfg.numrandomization    = 500;      % number of draws from the permutation distribution
cfg.frequency           = [5 30];
cfg.latency             = [-1.5 1.5];

[statTFUpvsDown] = ft_freqstatistics(cfg,  allSub_negPeak.up{:}, allSub_negPeak.down{:});
[statTFUpvsRest] = ft_freqstatistics(cfg,  allSub_negPeak.up{:}, allSub_negPeak.rest{:});
[statTFDownvsRest] = ft_freqstatistics(cfg,  allSub_negPeak.down{:}, allSub_negPeak.rest{:});

save([initPath.Exp '\data\group\eeg_stats\tfStat_UpVsDown.mat'],statTFUpvsDown) % for EEG12
save([initPath.Exp '\data\group\eeg_stats\tfStat_UpVsRest.mat'],statTFUpvsRest) % for EEG12
save([initPath.Exp '\data\group\eeg_stats\tfStat_DownVsRest.mat'],statTFDownvsRest) % for EEG12


% TF power for correlation (IRM)
power = [];

for idx_sub = 1 : length(allSub_negPeak.up)

    cfg = [];
    cfg.avgoverchan         = 'yes';
    cfg.latency             = [0.25 0.4];
    cfg.frequency           = [12 17];
    cfg.avgovertime         = 'yes';
    cfg.avgoverfreq         = 'yes';

    tmp = ft_selectdata(cfg,allSub_negPeak.up{idx_sub});

    power(idx_sub,1) = mean(tmp.powspctrm);
    
    tmp = ft_selectdata(cfg,allSub_negPeak.down{idx_sub});

    power(idx_sub,2) = mean(tmp.powspctrm);

  
end

csvwrite([initPath.Exp '\data\group\powerSigmaContrastUpvsDown.csv'],power)%=> to write getEEG() function



%% Effect sizes

% ERP
% Up vs Down
% + cluster
cfg = [];
cfg.channel = 'all';
cfg.latency = [0.32 0.64];
cfg.avgoverchan = 'yes';  
cfg.avgovertime = 'yes';  
up = ft_selectdata(cfg, grdAvgERP_negPeak.up);
down  = ft_selectdata(cfg, grdAvgERP_negPeak.down);

cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design(1,1:2*(length(listSubEEG)))  = [ones(1,(length(listSubEEG))) 2*ones(1,(length(listSubEEG)))];
cfg.design(2,1:2*(length(listSubEEG)))  = [1:(length(listSubEEG)) 1:(length(listSubEEG))];
effect_size.erp.UpvsDown.posCluster = ft_timelockstatistics(cfg, up, down);

% - cluster
cfg = [];
cfg.channel = 'all';
cfg.latency = [0.76 1.04];
cfg.avgoverchan = 'yes';  
cfg.avgovertime = 'yes';  
up = ft_selectdata(cfg, grdAvgERP_negPeak.up);
down  = ft_selectdata(cfg, grdAvgERP_negPeak.down);

cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design(1:2*(length(listSubEEG)),1)  = [ones(1,(length(listSubEEG))) 2*ones(1,(length(listSubEEG)))];
cfg.design(1:2*(length(listSubEEG)),2)  = [1:(length(listSubEEG)) 1:(length(listSubEEG))];
effect_size.erp.UpvsDown.negCluster = ft_timelockstatistics(cfg, up, down);

% Up vs Rest
% + cluster
cfg = [];
cfg.channel = 'all';
cfg.latency = [-0.61 0.62];
cfg.avgoverchan = 'yes';  
cfg.avgovertime = 'yes';  
up = ft_selectdata(cfg, grdAvgERP_negPeak.up);
rest  = ft_selectdata(cfg, grdAvgERP_negPeak.rest);

cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design(1,1:2*(length(listSubEEG)))  = [ones(1,(length(listSubEEG))) 2*ones(1,(length(listSubEEG)))];
cfg.design(2,1:2*(length(listSubEEG)))  = [1:(length(listSubEEG)) 1:(length(listSubEEG))];
effect_size.erp.UpvsRest.posCluster = ft_timelockstatistics(cfg, up, rest);

% - cluster
cfg = [];
cfg.channel = 'all';
cfg.latency = [0.72 1.14];
cfg.avgoverchan = 'yes';  
cfg.avgovertime = 'yes';  
up = ft_selectdata(cfg, grdAvgERP_negPeak.up);
rest  = ft_selectdata(cfg, grdAvgERP_negPeak.rest);

cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design(1,1:2*(length(listSubEEG)))  = [ones(1,(length(listSubEEG))) 2*ones(1,(length(listSubEEG)))];
cfg.design(2,1:2*(length(listSubEEG)))  = [1:(length(listSubEEG)) 1:(length(listSubEEG))];
effect_size.erp.UpvsRest.negCluster = ft_timelockstatistics(cfg, up, rest);

% Down vs Rest

% - cluster 2
cfg = [];
cfg.channel = 'all';
cfg.latency = [1.18 1.5];
cfg.avgoverchan = 'yes';  
cfg.avgovertime = 'yes';  
down = ft_selectdata(cfg, grdAvgERP_negPeak.down);
rest  = ft_selectdata(cfg, grdAvgERP_negPeak.rest);

cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design(1,1:2*(length(listSubEEG)))  = [ones(1,(length(listSubEEG))) 2*ones(1,(length(listSubEEG)))];
cfg.design(2,1:2*(length(listSubEEG)))  = [1:(length(listSubEEG)) 1:(length(listSubEEG))];
effect_size.erp.DownvsRest.negCluster2 = ft_timelockstatistics(cfg, down, rest);


% + cluster
cfg = [];
cfg.channel = 'all';
cfg.latency = [-0.64 0.42];
cfg.avgoverchan = 'yes';  
cfg.avgovertime = 'yes';  
down = ft_selectdata(cfg, grdAvgERP_negPeak.down);
rest  = ft_selectdata(cfg, grdAvgERP_negPeak.rest);

cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.ivar = 1;
cfg.uvar = 2;
cfg.design(1,1:2*(length(listSubEEG)))  = [ones(1,(length(listSubEEG))) 2*ones(1,(length(listSubEEG)))];
cfg.design(2,1:2*(length(listSubEEG)))  = [1:(length(listSubEEG)) 1:(length(listSubEEG))];
effect_size.erp.DownvsRest.posCluster = ft_timelockstatistics(cfg, down, rest);

% TF
% Up vs Down
[x,y]= find( squeeze(sum(statTFUpvsDown.mask)));
min(statTFUpvsDown.freq(x))
max(statTFUpvsDown.freq(x))
min(statTFUpvsDown.time(y))
max(statTFUpvsDown.time(y))
cfg = [];
cfg.channel = {'Fz' 'Cz' 'Pz' 'C3' 'C4'};
cfg.latency = [min(statTFUpvsDown.time(y)) max(statTFUpvsDown.time(y))];
cfg.frequency = [min(statTFUpvsDown.freq(x)) max(statTFUpvsDown.freq(x))];
cfg.avgoverchan = 'yes';   % this "squeezes" the channel dimension out of the data
cfg.avgovertime = 'yes';  % this "squeezes" the time dimension out of the data
cfg.avgoverfreq= 'yes';  % this "squeezes" the time dimension out of the data
statTFUpvsDown.roiUp    = ft_selectdata(cfg, grdAvg_negPeak.up);
statTFUpvsDown.roiDown  = ft_selectdata(cfg, grdAvg_negPeak.down);
cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.design(1,1:2*length(listSubEEG))  = [ones(1,length(listSubEEG)) 2*ones(1,length(listSubEEG))];
cfg.design(2,1:2*length(listSubEEG))  = [1:length(listSubEEG) 1:length(listSubEEG)];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
effect_size.tf.UpvsDown.posCluster = ft_freqstatistics(cfg, statTFUpvsDown.roiUp, statTFUpvsDown.roiDown);


% Up vs Rest
% - cluster
[x,y]= find( squeeze(sum((statTFUpvsRest.negclusterslabelmat==1) & statTFUpvsRest.mask))==5);
min(statTFUpvsRest.freq(x))
max(statTFUpvsRest.freq(x))
min(statTFUpvsRest.time(y))
max(statTFUpvsRest.time(y))

cfg = [];
cfg.latency = [min(statTFUpvsRest.time(y)) max(statTFUpvsRest.time(y))];
cfg.frequency = [min(statTFUpvsRest.freq(x)) max(statTFUpvsRest.freq(x))];
cfg.avgoverchan = 'yes';   % this "squeezes" the channel dimension out of the data
cfg.avgovertime = 'yes';  % this "squeezes" the time dimension out of the data
cfg.avgoverfreq= 'yes';  % this "squeezes" the time dimension out of the data
statTFUpvsRest.roiUp  = ft_selectdata(cfg, grdAvg_negPeak.up);
statTFUpvsRest.roiRest    = ft_selectdata(cfg, grdAvg_negPeak.rest);
cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.design(1,1:2*length(listSubEEG))  = [ones(1,length(listSubEEG)) 2*ones(1,length(listSubEEG))];
cfg.design(2,1:2*length(listSubEEG))  = [1:length(listSubEEG) 1:length(listSubEEG)];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

effect_size.tf.UpvsRest.negCluster = ft_freqstatistics(cfg, statTFUpvsRest.roiUp, statTFUpvsRest.roiRest);


% Down vs rest
% + cluster
[x,y]= find( squeeze(sum((statTFDownvsRest.posclusterslabelmat==1) & statTFDownvsRest.mask))==5);
min(statTFDownvsRest.freq(x))
max(statTFDownvsRest.freq(x))
min(statTFDownvsRest.time(y))
max(statTFDownvsRest.time(y))

cfg = [];
cfg.latency = [min(statTFDownvsRest.time(y)) max(statTFDownvsRest.time(y))];
cfg.frequency = [min(statTFDownvsRest.freq(x)) max(statTFDownvsRest.freq(x))];
cfg.avgoverchan = 'yes';   % this "squeezes" the channel dimension out of the data
cfg.avgovertime = 'yes';  % this "squeezes" the time dimension out of the data
cfg.avgoverfreq= 'yes';  % this "squeezes" the time dimension out of the data
statTFDownvsRest.roiDown  = ft_selectdata(cfg, grdAvg_negPeak.down);
statTFDownvsRest.roiRest    = ft_selectdata(cfg, grdAvg_negPeak.rest);
cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.design(1,1:2*length(listSubEEG))  = [ones(1,length(listSubEEG)) 2*ones(1,length(listSubEEG))];
cfg.design(2,1:2*length(listSubEEG))  = [1:length(listSubEEG) 1:length(listSubEEG)];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

effect_size.tf.DownvsRest.posCluster = ft_freqstatistics(cfg, statTFDownvsRest.roiDown, statTFDownvsRest.roiRest);


% - cluster
[x,y]= find( squeeze(sum((statTFDownvsRest.negclusterslabelmat==1) & statTFDownvsRest.mask))==5);
min(statTFDownvsRest.freq(x))
max(statTFDownvsRest.freq(x))
min(statTFDownvsRest.time(y))
max(statTFDownvsRest.time(y))

cfg = [];
cfg.latency = [min(statTFDownvsRest.time(y)) max(statTFDownvsRest.time(y))];
cfg.frequency = [min(statTFDownvsRest.freq(x)) max(statTFDownvsRest.freq(x))];
cfg.avgoverchan = 'yes';   % this "squeezes" the channel dimension out of the data
cfg.avgovertime = 'yes';  % this "squeezes" the time dimension out of the data
cfg.avgoverfreq= 'yes';  % this "squeezes" the time dimension out of the data
statTFDownvsRest.roiDown  = ft_selectdata(cfg, grdAvg_negPeak.down);
statTFDownvsRest.roiRest    = ft_selectdata(cfg, grdAvg_negPeak.rest);
cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'cohensd'; % see FT_STATFUN_COHENSD
cfg.design(1,1:2*length(listSubEEG))  = [ones(1,length(listSubEEG)) 2*ones(1,length(listSubEEG))];
cfg.design(2,1:2*length(listSubEEG))  = [1:length(listSubEEG) 1:length(listSubEEG)];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

effect_size.tf.DownvsRest.negCluster = ft_freqstatistics(cfg, statTFDownvsRest.roiDown, statTFDownvsRest.roiRest);

