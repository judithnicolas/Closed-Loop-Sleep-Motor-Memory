%Fig S10


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
    
    
    trl = load ([initPath.Exp '\data\' sub '\exp\eeg\' sub '_trl_trough_EEG_withFP.csv']) % time of negPeak / 1 = down ; 2 = rest; 3 = up; / sleep stage
    
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
cfg.keepindividual = 'no'; % for cohen's d computation
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

csvwrite([initPath.Exp '\data\group\nbTrialsEvoked_withFP.CSV'],allSub_negPeak.nbTrl)



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
cfg.alpha               = 0.05; 
cfg.clusteralpha        = 0.01; 
% cfg.tfce_H              = 2;       % default setting
% cfg.tfce_E              = 0.5;     % default setting
cfg.numrandomization    = 500;     
cfg.latency             = [-1.5 1.5];


[staterpUpVsDown] = ft_timelockstatistics(cfg,  allSubERP_negPeak.up{:}, allSubERP_negPeak.down{:});
[staterpUpVsRest] = ft_timelockstatistics(cfg,  allSubERP_negPeak.up{:}, allSubERP_negPeak.rest{:});
[staterpDownVsRest] = ft_timelockstatistics(cfg,  allSubERP_negPeak.down{:}, allSubERP_negPeak.rest{:});
save([initPath.Exp '\data\group\eeg_stats\erpStat_UpVsDown_withFP.mat'],staterpUpVsDown) % EEG_21
save([initPath.Exp '\data\group\eeg_stats\erpStat_UpVsRest_withFP.mat'],staterpUpVsRest) % EEG_21
save([initPath.Exp '\data\group\eeg_stats\erpStat_DownVsRest_withFP.mat'],staterpDownVsRest) % EEG_21


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
cfg.clusteralpha        = 0.01;
cfg.numrandomization    = 500;      % number of draws from the permutation distribution
cfg.frequency           = [5 30];
cfg.latency             = [-1.5 1.5];

cfg.alpha               = 0.05; 
[statTFUpvsDown] = ft_freqstatistics(cfg,  allSub_negPeak.up{:}, allSub_negPeak.down{:});
cfg.alpha               = 0.025/3; 
[statTFUpvsRest] = ft_freqstatistics(cfg,  allSub_negPeak.up{:}, allSub_negPeak.rest{:});
[statTFDownvsRest] = ft_freqstatistics(cfg,  allSub_negPeak.down{:}, allSub_negPeak.rest{:});

save([initPath.Exp '\data\group\eeg_stats\tfStat_UpVsDown_withFP.mat'],statTFUpvsDown) % EEG_21
save([initPath.Exp '\data\group\eeg_stats\tfStat_UpVsRest_withFP.mat'],statTFUpvsRest) % EEG_21
save([initPath.Exp '\data\group\eeg_stats\tfStat_DownVsRest_withFP.mat'],statTFDownvsRest) % EEG_21
