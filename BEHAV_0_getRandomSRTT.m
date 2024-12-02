%% get PVT from log file
% Judith Nicolas
% Created 2020 at KU Leuven

function BEHAV_0_getRandomSRTT(listSub,dirOutput,dirInput,nbSession,keyPresses )

    file= fopen([dirOutput 'randomSRTT.csv'],'w');
    fprintf(file,'%s; %s; %s; %s; %s; %s; %s; %s; %s\n', 'Sub' , 'Session','Block','Cue','timeCue','Rep','timeRep','RT','Acc');

    clc
    for idx_sub = 1 : length(listSub)
        sub = listSub{idx_sub};
        fprintf('%s\n',sub)

        behavFile = [dirInput sub '\behav\MSL_openLoop_' sub '.txt'];

        fid=fopen(behavFile);
        aline = fread(fid, 'char=>char');          % returns a single long string
        fclose(fid);

        aline(aline==uint8(sprintf('\r'))) = [];        % remove cariage return
        aline = tokenize(aline, uint8(sprintf('\n')));  % split on newline

        scanLines=find(~cellfun(@isempty,regexp(aline, '.*StartScan.*')));
        aline(scanLines)=[];        
        restLines=find(~cellfun(@isempty,regexp(aline, '.*rest.*')));
        aline(restLines)=[];        
        practiceLines=find(~cellfun(@isempty,regexp(aline, '.*Practice.*')));
        aline(practiceLines)=[];
        lines= find(~cellfun(@isempty,regexp(aline, '.*RC-.*')));
        if length(lines)/2> nbSession
            warning(['check logfile because more than 2 random SRTT have been launched\n Particpant to check: ' sub])
            %depends on design, suggestion is to have a raw log file and a
            %corrected logfile formated as if everything went smotth during
            %data acquisition 
        end
        
        lines= reshape(lines,2,length(lines)/2);
       
        nbBlock = (lines(2,1)-lines(1,1)-1)/2/keyPresses;

        for idx_sess= 1:size(lines,2)
            cue = aline(lines(1,idx_sess)+1:2:lines(2,idx_sess)-1);
            rep = aline(lines(1,idx_sess)+2:2:lines(2,idx_sess)-1);
            counter=1;

            
            sessionLine = strsplit(aline{lines(1,idx_sess)},'-');
            if strcmp(sessionLine{2}(1:2),'PR')
                session = 'preNight';
            elseif strcmp(sessionLine{2}(1:2),'PO')
                session = 'postNight';
            end

            
            for idx_block = 1 : nbBlock
               
                for idx_key = 1 : keyPresses
                    t0=cue{counter};t0=strsplit(t0,' ');
                    t1=rep{counter};t1=strsplit(t1,' ');
                    
  
                    if strcmp(t0{2},t1{2})
                        corr=1;
                    else
                        corr=0;
                    end
                    
                    fprintf(file,'%s; %s; %i; %s; %4.4f; %s; %4.4f; %4.4f; %i\n', ...
                        sub, session,idx_block,t0{2},str2num(t0{end}),t1{2},str2num(t1{end}),(str2num(t1{end})-str2num(t0{end}))*1000,corr);
                    counter=counter+1;
                end
            end
        end
    end
    fclose(file);
end