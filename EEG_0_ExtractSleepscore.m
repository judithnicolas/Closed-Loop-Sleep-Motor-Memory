%% TMR accuracy
% Judith Nicolas
% Created 2021 at KU Leuven
% Inspired from Genevieve Albouy


%% Incorporate score
pathIn = initPath.Exp;

[listSub,listSubEEG,listSubBehav] = getFullDatasets;

fileSleepDur = fopen([pathIn 'data\group\sleepDuration.csv'],'w');
fprintf(fileSleepDur,'%s; %s; %s; %s\n', 'Sub' , 'Stage','DV','Percentage');

fileSleepLat = fopen([pathIn 'data\group\sleepLatency.csv'],'w');
fprintf(fileSleepLat,'%s; %s; %s; %s; %s\n', 'Sub' , 'S1','S2', 'S3', 'REM');

fileSleepArousal = fopen([pathIn 'data\group\sleepArousal.csv'],'w');
fprintf(fileSleepArousal,'%s; %s; %s; %s\n', 'Sub' , 'Nb','perHour','Duration');

fileStimEff= fopen([pathIn 'data\group\stimEfficiency.csv'],'w');
fprintf(fileStimEff,'%s; %s; %s; %s; %s\n', 'Sub' , 'type','Stage','DV','Percentage'); 


for idx_sub = 1 : length(listSubEEG)
    
    sub = listSubEEG{idx_sub};
    fasstFile     = dir([pathIn '\data\' sub '\exp\eeg\' sub '.mat']); % obtained by opening the eeg file in FASST
    load([fasstFile.folder '/' fasstFile(1).name])
    eventfile     = dir([pathIn '\data\' sub '\exp\eeg\*.vmrk']);
    event = ft_read_event([eventfile.folder '/' eventfile(1).name]);


    scoreFile     = dir([pathIn '\data\' sub '\exp\score\*_ExtraitEvts.csv']); % Scores (_ExtraitEvts.csv) available at https://publicneuro.eu/catalogue.html
    score = tdfread(scoreFile(1).name,';');
    scoredEpochs=[];
    arousal=[];
    counterStade = 1;
    counterArousal = 1;
    for idx = 1 : length(score.Groupe)

        if strcmp(score.Groupe(idx,1),'S')
            scoredEpochs (counterStade,1) = score.Point_0x28Sec0xC9ch0x29(idx);
            scoredEpochs (counterStade,2) = score.Point_0x28Sec0xC9ch0x29(idx)+score.Dur0xE9e_0x28Points0x29(idx)-1;
            scoredEpochs (counterStade,3) = score.Stade(idx);
            counterStade = counterStade +1;
        elseif strcmp(score.Groupe(idx,1),'M')
            arousal (counterArousal,1) = score.Point_0x28Sec0xC9ch0x29(idx);
            arousal (counterArousal,2) = score.Dur0xE9e_0x28Points0x29(idx);
            counterArousal = counterArousal +1;
        end

    end

    scoredEpochs(1,1) = 1;
    scoredEpochs(end,2) = D.Nsamples-1;


    D.other.CRC.score=[];
    D.other.CRC.score{1}= scoredEpochs(:,3)' ;
    D.other.CRC.score(2,1)= {'Sonia'};
    D.other.CRC.score(3,1)= {30};
    D.other.CRC.score{4,1} = [0.001 D.Nsamples/D.Fsample];
    D.other.CRC.score{5,1} = []; D.other.CRC.score{7,1} = []; D.other.CRC.score{8,1} = []; 
    if isempty(arousal)
        D.other.CRC.score{6,1}=[];
    else
        D.other.CRC.score{6,1}=[arousal(:,1)/D.Fsample arousal(:,1)/D.Fsample+arousal(:,2)/D.Fsample];
    end


    selStims  = strcmp({event.value}, {'S 29'});
    scoredStims= [];
    for idx_event = 1 : length(selStims )

        if selStims(idx_event)==1
            f1 = (event(idx_event).sample>scoredEpochs(:,1));
            f2 = (event(idx_event).sample < scoredEpochs(:,2));
            epoch_stim = find(f1.*f2);

            scoredStims = vertcat(scoredStims,[event(idx_event).sample scoredEpochs(epoch_stim,3)]);
        end
    end
    save([pathIn '\data\' sub '\exp\eeg\' sub '_scoredStims.mat'],'scoredStims')
    save([pathIn '\data\' sub '\exp\eeg\' sub '_scoredEpochs.mat'],'scoredEpochs')


    
    event = ft_read_event([eventfile.folder '/' eventfile.name]);
    %     compute sleep metrics
    %     Determine ODL&CDL
    
    doorLightMkr = [];
   
    if sum(strcmp({event.value}, {'cdl'}))  >0
        selCDL  = strcmp({event.value}, {'cdl'});
        doorLightMkr(1) = event(selCDL).sample;
    else
        doorLightMkr(1) = 0;
    end
 
    if sum(strcmp({event.value}, {'odl'}))>0
        selODL  = strcmp({event.value}, {'odl'});
        doorLightMkr(2) = event(selODL  ).sample;
    elseif sum(strcmp({event.value}, {'ODL'}))>0
        selODL    = strcmp({event.value}, {'ODL'});
        doorLightMkr(2) = event(selODL  ).sample;
    else
        doorLightMkr(2) = D.Nsamples/D.Fsample*1000;
        
    end
   

    scoredStims  = scoredStims(scoredStims (:,1)>doorLightMkr(1) & scoredStims (:,1)<doorLightMkr(2),:)  ;
    scoredStims(:,3) = scoredStims(:,2);
    
    WakeCues            = length(find(scoredStims(:,2) == 0));
    S1Cues              = length(find(scoredStims(:,2) == 1));
    S2Cues              = length(find(scoredStims(:,2) == 2));
    S3Cues              = length(find(scoredStims(:,2) == 3));
    REMCues             = length(find(scoredStims(:,2) == 5));
    accurateCues        = S2Cues+S3Cues;
    totCues             = length(scoredStims);
    accuracy            = (accurateCues/totCues)*100;
    fprintf('%s\t%.2f\n',sub,accuracy)
    
    WakeUp            = length(find(scoredStims(:,2) == 0 & scoredStims(:,3) == 1));
    S1Up              = length(find(scoredStims(:,2) == 1 & scoredStims(:,3) == 1));
    S2Up              = length(find(scoredStims(:,2) == 2 & scoredStims(:,3) == 1));
    S3Up              = length(find(scoredStims(:,2) == 3 & scoredStims(:,3) == 1));
    REMUp             = length(find(scoredStims(:,2) == 5 & scoredStims(:,3) == 1));
    TotalUp           = length(find(scoredStims(:,3) == 1));
     
    WakeDown          = length(find(scoredStims(:,2) == 0 & scoredStims(:,3) == -1));
    S1Down            = length(find(scoredStims(:,2) == 1 & scoredStims(:,3) == -1));
    S2Down            = length(find(scoredStims(:,2) == 2 & scoredStims(:,3) == -1));
    S3Down            = length(find(scoredStims(:,2) == 3 & scoredStims(:,3) == -1));
    REMDown           = length(find(scoredStims(:,2) == 5 & scoredStims(:,3) == -1));
    TotalDown         = length(find(scoredStims(:,3) == -1));
    

    %% Extract sleep characteristics
    % Create variables
    TRS = []; % Time allowed to sleep (Between CDL (close door light) et ODL (Open door light))
    TPS = []; % Time of the sleeping period
    TST = []; % Total Sleep Period
    LatS1 = []; % S1 Latency
    LatS2 = []; % S2 Latency
    LatREM = []; % REM Latency
    W = []; % Time awake
    S1 = []; % Time in S1
    S2 = []; % Time in S2
    S3 = []; % Time in S3
    REM = []; % Time in REM
    MT = []; % Time in MT
    SEff = []; % Sleep Efficiency
    S1Eff = []; % S1 Efficiency
    S2Eff = [];% S2 Efficiency
    S3Eff = [];% S3 Efficiency
    S4Eff = [];% S4 Efficiency
    REMEff = [];% REM Efficiency
    nbar = []; % Number of Arousal
    Arhour = []; % Number of Arousal per hour
    Ardur = []; % Mean Arousal duration
    
    
    Winsize = D.other.CRC.score{3,size(D.other.CRC.score,2)};
    
    %Time allowed to sleep
    TRS = [doorLightMkr(2) - doorLightMkr(1)]/1000/60;

    %Invalidation
    adapted = 1:length(D.other.CRC.score{1,size(D.other.CRC.score,2)});
    nottobescored = find(adapted < doorLightMkr(1)/Winsize | adapted > doorLightMkr(2)/Winsize);
    iisc = D.other.CRC.score{1,size(D.other.CRC.score,2)};
    iisc(nottobescored)=-1;
    
    Zero  = find(scoredEpochs(:,3) == 0 & scoredEpochs (:,1)>doorLightMkr(1) & scoredEpochs (:,2)<doorLightMkr(2));
    One   = find(scoredEpochs(:,3) == 1 & scoredEpochs (:,1)>doorLightMkr(1) & scoredEpochs (:,2)<doorLightMkr(2));
    Two   = find(scoredEpochs(:,3) == 2 & scoredEpochs (:,1)>doorLightMkr(1) & scoredEpochs (:,2)<doorLightMkr(2));
    Three = find(scoredEpochs(:,3) == 3 & scoredEpochs (:,1)>doorLightMkr(1) & scoredEpochs (:,2)<doorLightMkr(2));
    Five  = find(scoredEpochs(:,3) == 5 & scoredEpochs (:,1)>doorLightMkr(1) & scoredEpochs (:,2)<doorLightMkr(2));
    
    
    %Total Sleep Time
    TST = [length([Two ;Three ;Five])*Winsize]/60;
    
    %Lantency St1
    if length(One) == 0
        LatS1=nan;
    else
        LatS1 = [One(1)*Winsize - doorLightMkr(1)/1000]/60;
    end
    
    %Lantency St2
    if length(Two) == 0
        LatS2=nan;
    else
        LatS2 = [Two(1)*Winsize - doorLightMkr(1)/1000]/60;
    end
    
    %Lantency St3
    if length(Three) == 0
        LatS3=nan;
    else
        LatS3 = [Three(1)*Winsize - doorLightMkr(1)/1000]/60;
    end
    
    %Lantency REM
    if length(Five) == 0
        LatREM = nan;
    else
        LatREM = [Five(1)*Winsize - doorLightMkr(1)/1000]/60;
    end
    
    %Min & Percentage of W
    W = [length([Zero])*Winsize]/60;
    
    %Min & Percentage of S1
    S1 = [length([One])*Winsize]/60;
    
    %Min & Percentage of S2
    S2 = [length([Two])*Winsize]/60;
    
    %Min & Percentage of S3
    S3 =[ length([Three])*Winsize]/60;
    
    %Min & Percentage of REM
    REM = [length([Five])*Winsize]/60;
    
    
    % Extract micro-arousals
    if isempty(D.other.CRC.score{6,1})
        nBar = 0;Arhour= 0; Ardur=0;
    else
        nBar = size(D.other.CRC.score{6,1},1);
        Arhour = nBar/(TST/60); %TST in hour
        Ardur = sum(D.other.CRC.score{6,1}(:,2)-D.other.CRC.score{6,1}(:,1)); %in seconds
    end
    
    fprintf(fileSleepArousal,'%s; %i; %5.3f; %5.3f\n',sub, nBar,Arhour,Ardur);
    
    
    fprintf(fileSleepDur,'%s; %s; %5.3f; %3.2f\n', sub , 'opportunity',TRS,TRS/TRS*100);
    fprintf(fileSleepDur,'%s; %s; %5.3f; %3.2f\n', sub , 'sleep',TST,TST/TRS*100);
    fprintf(fileSleepDur,'%s; %s; %5.3f; %3.2f\n', sub , 'wake',W,W/TRS*100);
    fprintf(fileSleepDur,'%s; %s; %5.3f; %3.2f\n', sub , 'S1',S1,S1/TRS*100);
    fprintf(fileSleepDur,'%s; %s; %5.3f; %3.2f\n', sub , 'S2',S2,S2/TRS*100);
    fprintf(fileSleepDur,'%s; %s; %5.3f; %3.2f\n', sub , 'S3',S3,S3/TRS*100);
    fprintf(fileSleepDur,'%s; %s; %5.3f; %3.2f\n', sub , 'REM',REM,REM/TRS*100);
    fprintf(fileSleepLat,'%s; %5.3f; %5.3f; %5.3f; %5.3f\n', sub , LatS1,LatS2, LatS3, LatREM);
    
    
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'all','total',totCues,totCues/totCues*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'all','wake',WakeCues,WakeCues/totCues*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'all','S1',S1Cues,S1Cues/totCues*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'all','S2',S2Cues,S2Cues/totCues*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'all','S3',S3Cues,S3Cues/totCues*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'all','REM',REMCues,REMCues/totCues*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'all','NREM',S2Cues+S3Cues,(S2Cues+S3Cues)/totCues*100);
    
    
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'up','total',TotalUp,TotalUp/totCues*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'up','wake',WakeUp,WakeUp/TotalUp*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'up','S1',S1Up,S1Up/TotalUp*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'up','S2',S2Up,S2Up/TotalUp*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'up','S3',S3Up,S3Up/TotalUp*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'up','REM',REMUp,REMUp/TotalUp*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'up','NREM',S2Up+S3Up,(S2Up+S3Up)/TotalUp*100);
    
    
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'down','total',TotalDown,TotalDown/totCues*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'down','wake',WakeDown,WakeDown/TotalDown*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'down','S1',S1Down,S1Down/TotalDown*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'down','S2',S2Down,S2Down/TotalDown*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'down','S3',S3Down,S3Down/TotalDown*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'down','REM',REMDown,REMDown/TotalDown*100);
    fprintf(fileStimEff,'%s; %s; %s; %i; %3.2f\n', sub , 'down','NREM',S2Down+S3Down,(S2Down+S3Down)/TotalDown*100);
    
end


fclose(fileSleepDur);
fclose(fileStimEff);
fclose(fileSleepArousal);
fclose(fileSleepLat);
% %%%END OF SCRIPT
