%% get PVT from log file
% Judith Nicolas
% Created 2020 at KU Leuven

function BEHAV_0_getSequentialSRTT(listSub,dirInput,dirOutput)

run parameters.m


filetot = fopen([dirOutput '\sequentialSRTT_tmp.csv'],'w');
fprintf(filetot,'%s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s; %s\n', 'Sub' , 'Session','Block','Sequence','Condition', 'Sound', 'Order','OrdinalPos','Repetition','Cue','timeCue','Rep','timeRep','RT','StartScan');

clc
for idx_sub = 1 : length(listSub)
    sub = listSub{idx_sub};
    
    behavFile = [dirInput sub '\behav\MSL_openLoop_' sub '.txt'];
    load( [dirInput sub '\behav\MSL_openLoop_' sub '.mat']);

    fid=fopen(behavFile);
    aline = fread(fid, 'char=>char');          % returns a single long string
    fclose(fid);
    
    aline(aline==uint8(sprintf('\r'))) = [];        % remove cariage return
    aline = tokenize(aline, uint8(sprintf('\n')));  % split on newline
    scanLines_idx=find(~cellfun(@isempty,regexp(aline, '.*StartScan.*')));
    
    scanTimes = aline(scanLines_idx);
    
    if length(scanTimes)>5
        warning('Check text file for multiple sessions')
    end
    startScan=[];
    for idx_line = scanLines_idx
        taskMSL = strsplit(aline{idx_line-1});
        taskMSL = taskMSL{end-1};
        time = strsplit(aline{idx_line});time = str2num(time{end});
        
        if strcmp(taskMSL,'SL-TRAININGPRE_NIGHT')
            startScan(1) = time;
        
        elseif strcmp(taskMSL,'SL-TEST-PRE_NIGHT')
            startScan(2) = time;
        elseif strcmp(taskMSL,'SL-TRAINING-POST_NIGHT')
            startScan(3) = time;
        end
    end
    
    
    aline(scanLines_idx)=[];    
    auditoryCueLines=find(~cellfun(@isempty,regexp(aline, '.*autory.*')));
    aline(auditoryCueLines)=[];    
    restLines=find(~cellfun(@isempty,regexp(aline, '.*rest.*')));
    aline(restLines)=[];
    practiceLines=find(~cellfun(@isempty,regexp(aline, '.*Practice.*')));
    aline(practiceLines)=[];
    lines= find(~cellfun(@isempty,regexp(aline, '.*SL-.*')));
    lines= reshape(lines,2,length(lines)/2);
    
    lines(1,:)=lines(1,:)+1;
    lines(2,:)=lines(2,:)-1;
    
    cue             = [];
    rep             = [];
    counterBlock    = 1;
    block           = 1;
    counterOrdinal  = 1;
    repetition      = 1;
    for idx_sess = 1 : size(lines,2)
        
        allLines = aline(lines(1,idx_sess):lines(2,idx_sess))';
        
        for idx_line = 1 : length(allLines)
            
            tmpCue         = strsplit(allLines{idx_line});
            
            if strcmp(tmpCue{3},'cue') 
                cueKey      = str2num(tmpCue{4});
                cueTime     = str2num(tmpCue{end});

                tmpResp     = strsplit(allLines{idx_line+1});
                respKey     = str2num(tmpResp{4});
                respTime    = str2num(tmpResp{end});


                if strcmp(sub,'CL_04') && strcmp(seq,'B')  && counterOrdinal == 4
                else
                    if cueKey == param.sequence(1,counterOrdinal)
                        seq     = 'A';
                        cond    = param.Seq1.condition ;
                        sound   = param.Seq1.sound;
                        order   = find(param.seqOrder == 1);
                    elseif cueKey == param.sequence(2,counterOrdinal)
                        seq     = 'B';
                        cond    = param.Seq2.condition ;
                        sound   = param.Seq2.sound;
                        order   = find(param.seqOrder == 2);
                    elseif cueKey == param.sequence(3,counterOrdinal)
                        seq     = 'C';
                        cond    = param.Seq3.condition ;
                        sound   = param.Seq3.sound;
                        order   = find(param.seqOrder == 3);

                    end
                end
                if ismember(counterBlock,[1:5 21:25 41:45])
                    repetition = 1;
                elseif ismember(counterBlock,[6:10 26:30 46:50])
                    repetition = 2;
                elseif ismember(counterBlock,[11:15 31:35 51:55])
                    repetition = 3;
                elseif ismember(counterBlock,[16:20 36:40 56:60])
                    repetition = 4;
                end

                fprintf(filetot,'%s; %i; %i; %s; %s; %i; %i; %i; %i; %i; %4.4f; %i; %4.4f; %4.4f; %4.4f\n', ...
                    sub, idx_sess,block,seq,cond,sound,order,counterOrdinal,repetition,cueKey,cueTime,respKey,respTime,respTime-cueTime,startScan(idx_sess));

                if counterOrdinal < size(param.sequence,2)
                    counterOrdinal  = counterOrdinal +1;
                else
                    counterOrdinal  = 1;
                end

                if counterBlock < param.nbKeys*size(param.sequence,1)
                    counterBlock    = counterBlock +1;
                else
                    counterBlock    = 1;
                    block           = block +1;
                end

            end
            
        end
        
    end
end
fclose(filetot);
