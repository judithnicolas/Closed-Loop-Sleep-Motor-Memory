%% Load data & baseline
clc
[listSub,listSubEEG,listSubBehav] = getFullDatasets;

allSubERP = {};
allSub = {};

for idx_sub = 1 : length(listSubEEG)

    sub = listSubEEG{idx_sub};
    load ([initPath.Exp '\data\' sub '\exp\eeg\' sub '_preprocessed_continuous.mat'])
    trl = load ([initPath.Exp '\data\' sub '\exp\eeg\' sub '_trl_sham_onlyTP.csv']);
    disp (['loading ' sub ' dataset'])

    cfg=[];
    cfg.trl    = [trl(:,1)-2.5*data.fsample trl(:,1)+2*data.fsample repmat(-2.5*data.fsample,length(trl(:,1)),1)];
    data_epoched = ft_redefinetrial(cfg, data);

    cfg = [];
    cfg.channel        = [1:6];
    cfg.resamplefs      = 100;
    cfg.demean          = 'yes';
    data_epoched = ft_resampledata(cfg, data_epoched);

    cfg = [];
    cfg.channel        = [1:6];
    cfg.baselinewindow = [-2.5 -2];
    cfg.demean         = 'yes';
    data_epoched = ft_preprocessing(cfg, data_epoched);


    data_epoched.trialinfo = trl;

%%%%%%%%%%%%%%%%%%%%%%%%%
    % ERP
    % NREM up vs down and all

    % true up
    cfg = [];
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,2)== 3);
    timelock_data= ft_timelockanalysis(cfg, data_epoched);

    allSubERP.nbTrl(idx_sub,1) =  length(cfg.trials);
    timelock_data.cfg=[];
    allSubERP.trueUp{idx_sub}= timelock_data;

    %true down
    cfg = [];
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,2)== 1 );
    timelock_data= ft_timelockanalysis(cfg, data_epoched);
    allSubERP.nbTrl(idx_sub,2) =  length(cfg.trials);
    timelock_data.cfg=[];
    allSubERP.trueDown{idx_sub}= timelock_data;

    % sham up
    cfg = [];
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,2)== -3);
    timelock_data= ft_timelockanalysis(cfg, data_epoched);

    allSubERP.nbTrl(idx_sub,3) =  length(cfg.trials);
    timelock_data.cfg=[];
    allSubERP.shamUp{idx_sub}= timelock_data;

    % sham down
    cfg = [];
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,2)== -1 );
    timelock_data= ft_timelockanalysis(cfg, data_epoched);
    allSubERP.nbTrl(idx_sub,4) =  length(cfg.trials);
    timelock_data.cfg=[];
    allSubERP.shamDown{idx_sub}= timelock_data;


%%%%%%%%%%%%%%%%%%%%%%%%%
    % TF

        % NREM up vs down and all

    foi = 5:0.5:30; % 0.1 Hz steps
    toi = -2:0.01:3; % 0.1 s steps

    % true up
    cfg = [];
    cfg.method     = 'mtmconvol';
    cfg.pad        = 'nextpow2';
    cfg.taper      = 'hanning';
    cfg.foi        = foi;
    cfg.toi        = toi;
    cfg.t_ftimwin  = 5./cfg.foi;
    cfg.tapsmofrq  = 0.4 *cfg.foi;   
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,2)==3);
    timelock_data  = ft_freqanalysis(cfg, data_epoched);
    allSub.nbTrl(idx_sub,1) =  length(cfg.trials);

    cfg = [];
    cfg.baseline     =    [-2.5 -2];
    cfg.baselinetype = 'relchange';
    timelock_data    = ft_freqbaseline(cfg, timelock_data);

    timelock_data.cfg        =[];
    timelock_data.freq       = foi;
    allSub.trueUp{idx_sub}= timelock_data;

    % true down
    cfg = [];
    cfg.method     = 'mtmconvol';
    cfg.pad        = 'nextpow2';
    cfg.taper      = 'hanning';
    cfg.foi        = foi; % 0.1 Hz steps
    cfg.toi        = toi; % 0.1 s steps
    cfg.t_ftimwin  = 5./cfg.foi;
    cfg.tapsmofrq  = 0.4 *cfg.foi;    
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,2)==1 );
    timelock_data  = ft_freqanalysis(cfg, data_epoched);

    allSub.nbTrl(idx_sub,2) =  length(cfg.trials);

    cfg = [];
    cfg.baseline     =    [-2.5 -2];
    cfg.baselinetype = 'relchange';
    timelock_data    = ft_freqbaseline(cfg, timelock_data);

    timelock_data.cfg        =[];
    timelock_data.freq       = foi;
    allSub.trueDown{idx_sub}= timelock_data;


    % sham up
    cfg = [];
    cfg.method     = 'mtmconvol';
    cfg.pad        = 'nextpow2';
    cfg.taper      = 'hanning';
    cfg.foi        = foi;
    cfg.toi        = toi;
    cfg.t_ftimwin  = 5./cfg.foi;
    cfg.tapsmofrq  = 0.4 *cfg.foi;   
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,2)==-3);
    timelock_data  = ft_freqanalysis(cfg, data_epoched);
    allSub.nbTrl(idx_sub,3) =  length(cfg.trials);

    cfg = [];
    cfg.baseline     =    [-2.5 -2];
    cfg.baselinetype = 'relchange';
    timelock_data    = ft_freqbaseline(cfg, timelock_data);

    timelock_data.cfg        =[];
    timelock_data.freq       = foi;
    allSub.shamUp{idx_sub}= timelock_data;

    % sham down
    cfg = [];
    cfg.method     = 'mtmconvol';
    cfg.pad        = 'nextpow2';
    cfg.taper      = 'hanning';
    cfg.foi        = foi; % 0.1 Hz steps
    cfg.toi        = toi; % 0.1 s steps
    cfg.t_ftimwin  = 5./cfg.foi;
    cfg.tapsmofrq  = 0.4 *cfg.foi;    
    cfg.keeptrials = 'no';
    cfg.trials     = find(data_epoched.trialinfo(:,2)==-1 );
    timelock_data  = ft_freqanalysis(cfg, data_epoched);

    allSub.nbTrl(idx_sub,4) =  length(cfg.trials);

    cfg = [];
    cfg.baseline     =    [-2.5 -2];
    cfg.baselinetype = 'relchange';
    timelock_data    = ft_freqbaseline(cfg, timelock_data);

    timelock_data.cfg        =[];
    timelock_data.freq       = foi;
    allSub.shamDown{idx_sub}= timelock_data;

    % subtration erp
    cfg = [];
    cfg.parameter = 'avg';
    cfg.operation = 'x1-x2';
    allSubERP.subtractedUp{idx_sub}= ft_math(cfg,allSubERP.trueUp{idx_sub},allSubERP.shamUp{idx_sub});
    allSubERP.subtractedDown{idx_sub}= ft_math(cfg,allSubERP.trueDown{idx_sub},allSubERP.shamDown{idx_sub});

    % subtration tf
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'x1-x2';
    allSub.subtractedUp{idx_sub}= ft_math(cfg,allSub.trueUp{idx_sub},allSub.shamUp{idx_sub});
    allSub.subtractedDown{idx_sub}= ft_math(cfg,allSub.trueDown{idx_sub},allSub.shamDown{idx_sub});
    allSub.subtractedUpvsDown{idx_sub}= ft_math(cfg,allSub.trueUp{idx_sub},allSub.trueDown{idx_sub});

end

%% Avg
cfg = [];
% cfg.keepindividual = 'yes'; % for cohen's d computation
grdAvgERP.trueUp = ft_timelockgrandaverage(cfg, allSubERP.trueUp{:});
grdAvgERP.trueDown = ft_timelockgrandaverage(cfg, allSubERP.trueDown{:});
grdAvgERP.shamUp= ft_timelockgrandaverage(cfg, allSubERP.shamUp{:});
grdAvgERP.shamDown = ft_timelockgrandaverage(cfg, allSubERP.shamDown{:});
grdAvgERP.substractedUp = ft_timelockgrandaverage(cfg, allSubERP.subtractedUp{:});
grdAvgERP.substractedDown= ft_timelockgrandaverage(cfg, allSubERP.subtractedDown{:});

grdAvg.trueUp = ft_freqgrandaverage(cfg, allSub.trueUp{:});
grdAvg.trueDown = ft_freqgrandaverage(cfg, allSub.trueDown{:});
grdAvg.shamUp = ft_freqgrandaverage(cfg, allSub.shamUp{:});
grdAvg.shamDown = ft_freqgrandaverage(cfg, allSub.shamDown{:});
grdAvg.subtractedUp = ft_freqgrandaverage(cfg, allSub.subtractedUp{:});
grdAvg.subtractedDown = ft_freqgrandaverage(cfg, allSub.subtractedDown{:});
grdAvg.subtractedUpvsDown = ft_freqgrandaverage(cfg, allSub.subtractedUpvsDown{:});

clearvars -except initPath grdAvg allSub allSubERP  grdAvgERP listSubEEG


%% stat erp
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
cfg.alpha               = 0.05/6; 
cfg.clusteralpha        = 0.01; 
% cfg.tfce_H              = 2;       % default setting
% cfg.tfce_E              = 0.5;     % default setting
cfg.numrandomization    = 500;     
cfg.latency             = [-1.5 1.5];


[staterpUp] = ft_timelockstatistics(cfg,  allSubERP.trueUp{:}, allSubERP.shamUp{:});
[staterpDown] = ft_timelockstatistics(cfg,  allSubERP.trueDown{:}, allSubERP.shamDown{:});
[staterpupvsDown] = ft_timelockstatistics(cfg,  allSubERP.trueUp{:}, allSubERP.trueDown{:});

save([initPath.Exp '\data\group\eeg_stats\erpStat_trueUpVsShamUp.mat'],staterpUp) % from EEG_11_troughLockedAnalysis.m
save([initPath.Exp '\data\group\eeg_stats\erpStat_trueDownVsShamDown.mat'],staterpDown) % from EEG_11_troughLockedAnalysis.m
save([initPath.Exp 'data\group\eeg_stats\erpStat_trueUpvsTrueDown.mat'],staterpupvsDown) % from EEG_11_troughLockedAnalysis.m


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
cfg.alpha               = 0.05/6;
cfg.clusteralpha        = 0.01;
cfg.numrandomization    = 500;      % number of draws from the permutation distribution
cfg.frequency           = [5 30];
cfg.latency             = [-1.5 1.5];

[statTFUp] = ft_freqstatistics(cfg,  allSub.trueUp{:}, allSub.shamUp{:});
[statTFDown] = ft_freqstatistics(cfg,  allSub.trueDown{:}, allSub.shamDown{:});
[statTFupvsDown] = ft_freqstatistics(cfg,  allSub.trueUp{:}, allSub.trueDown{:});

save([initPath.Exp '\data\group\eeg_stats\tfStat_trueUpVsShamUp.mat'],statTFUp) % from EEG_11_troughLockedAnalysis.m
save([initPath.Exp '\data\group\eeg_stats\tfStat_trueDownVsShamDown.mat'],statTFDown) % from EEG_11_troughLockedAnalysis.m
save([initPath.Exp 'data\group\eeg_stats\tfStat_trueUpvsTrueDown.mat'],statTFupvsDown) % from EEG_11_troughLockedAnalysis.m
