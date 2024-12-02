% Pre-proccess sleep EEG data based on scored epochs (from Extract_Sleepscore.m)

[listSub,listSubEEG,listSubBehav] = getFullDatasets;

for idx_sub = 1 : length(listSubEEG)
    
    sub = listSubEEG{idx_sub};
    fprintf('%s\n',sub)

    load([initPath.Exp '\data\' sub '\exp\eeg\' sub '_scoredEpochs.mat'])
    load([initPath.Exp '\data\' sub '\exp\eeg\' sub '.mat']) % from extract_SleepScore.m

    tmp = dir([initPath.Exp '\data\' sub '\exp\eeg\' sub '*.edf']);

    % re-ref
    cfg = [];
    cfg.trl = scoredEpochs(scoredEpochs(:,3)==2 |scoredEpochs(:,3)==3,: );
    cfg.dataset = [initPath.Exp '\data\' sub '\exp\eeg\' tmp.name]; 
    cfg.channel     = 'all';
    cfg.reref       = 'yes'; % not for CL_10/CL_13/CL_15/CL_22/CL_29 because of edge artifacts
    cfg.implicitref = 'REF1';        
    cfg.refchannel  = {'REF1', 'REF2'}; 
    cfg.continuous   = 'yes' ;
    data_reref= ft_preprocessing(cfg);
   
        
    % Filtering
    cfg = [];
    cfg.channel       = 1:6;
    cfg.hpfilter      = 'yes'; % highpass filter (default = 'no')
    cfg.hpfreq        = 0.1;   % highpass frequency in Hz 
    cfg.hpfiltord     = 4;   % filter order
    cfg.lpfilter      = 'yes'; % lowpass filter (default = 'no')
    cfg.lpfreq        = 30;    % lowpass  frequency in Hz
    cfg.continuous    = 'yes';
    data_EEG_filtered = ft_preprocessing(cfg, data_reref);
    
    cfg = [];
    cfg.channel       = 7:10;
    cfg.hpfilter      = 'yes'; % highpass filter (default = 'no')
    cfg.hpfreq        = 2;   % highpass frequency in Hz
    cfg.lpfilter      = 'yes'; % lowpass filter (default = 'no')
    cfg.lpfreq        = 15;    % lowpass  frequency in Hz
    cfg.continuous    = 'yes';
    data_EOG_filtered = ft_preprocessing(cfg, data_reref);
    
    cfg = [];
    cfg.channel       = 11:12;
    cfg.hpfilter      = 'yes'; % highpass filter (default = 'no')
    cfg.hpfreq        = 100;   % highpass frequency in Hz
    cfg.lpfilter      = 'yes'; % lowpass filter (default = 'no')
    cfg.lpfreq        = 150;    % lowpass  frequency in Hz
    cfg.continuous    = 'yes';
    data_EMG_filtered = ft_preprocessing(cfg, data_reref);
    
    cfg = [];
    data_filtered= ft_appenddata(cfg,data_EEG_filtered,data_EOG_filtered,data_EMG_filtered);
    
    cfg=[];
    cfg.artfctdef.reject = 'nan'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
    cfg.artfctdef.muscle.artifact = D.other.CRC.score{6,1}*D.Fsample;
    data_arousal_artifacts = ft_rejectartifact(cfg,data_filtered);
    
    
    % Artefact inspection Data browser
    cfg = [];
    cfg.viewmode   = 'vertical';
    cfg.ylim       = [-150 150];
    cfg.continuous = 'yes';
    cfg.blocksize  = 30;
    cfg.channel    = [1:3 5:6 11:12];
    artf=ft_databrowser(cfg,data_arousal_artifacts);
    
    cfg=[];
    cfg.artfctdef.reject = 'nan'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
    cfg.artfctdef.visual.artifact =  artf.artfctdef.visual.artifact;
    cfg.artfctdef.muscle.artifact = D.other.CRC.score{6,1}*D.Fsample;
    data_artifacts = ft_rejectartifact(cfg,data_filtered);
    
    save ([initPath.Exp '\data\' sub '\exp\eeg\' sub '_data_artifacts_new.mat'], 'data_artifacts')

        
end


%% ICA Computation
    

for idx_sub = 1 : length(listSubEEG)
    
    sub = listSubEEG{idx_sub};
    fprintf('%s\n',sub)

    load ([initPath.Exp '\data\' sub '\exp\eeg\' sub '_data_artifacts.mat'])

    cfg = [];
    ic_data = ft_componentanalysis(cfg,data);
    ic_data = ft_componentanalysis(cfg,data_artifacts);
    
    save ([initPath.Exp '\data\' sub '\exp\eeg\' sub '_ICA_data.mat'], 'ic_data')
    
    
end

%% ICA inspection

for idx_sub = 1 : length(listSubEEG) 
    
    sub = listSubEEG{idx_sub};

    fprintf('%s\n',sub)

    load ([initPath.Exp '\data\' sub '\exp\eeg\' sub '_ICA_data.mat'])
    load ([initPath.Exp '\data\' sub '\exp\eeg\' sub '_data_artifacts.mat'])

    cfg          = [];
    cfg.viewmode = 'component';
    cfg.layout =  [initPath.FieldTrip '\template\layout\EEG1020.lay'];
    ft_databrowser(cfg, ic_data)
    
    [components] = input('Value(s) for component ([x1 x2...] if more than one values): ');
    
    cfg = [];
    cfg.component = components; % to be removed component(s)
    data = ft_rejectcomponent(cfg, ic_data, data_artifacts);
    
    cfg=  [];
    cfg.channel = 1:6;
    data = ft_selectdata(cfg,data);
    
    save ([initPath.Exp '\data\' sub '\exp\eeg\' sub '_preprocessed_continuous.mat'], 'data')
    
        
    delete ([initPath.Exp '\data\' sub '\exp\eeg\' sub '_ICA_data.mat'])
    delete ([initPath.Exp '\data\' sub '\exp\eeg\' sub '_data_artifacts.mat'])

    
end
