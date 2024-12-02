%% get PVT from log file
% Judith Nicolas
% Created 2021 at KU Leuven

function BEHAV_0_getPVT(listSub,dirOutput,dirInput)

    file = fopen([dirOutput 'PVT.csv'],'w');
    fprintf(file,'%s; %s; %s; %s\n', 'Sub' ,'Session','Trial','RT');
    clc
    for idx_sub = 1 : length(listSub)
        sub = listSub{idx_sub};
        fprintf('%s\n',sub)
        
        if strcmp(sub,'CL_17')% lost
            fprintf(file,'%s; %s; %i; %s\n',sub, session,counter,'NA');
            for idx_sess= 1:2
                for idx_trl = 1:100
                    
                    if idx_sess== 1
                        session = 'preNight';
                    elseif idx_sess==2
                        session = 'postNight';
                    end
                
                fprintf(file,'%s; %s; %i; %s\n',sub, session,idx_trl,'NA');

                end
            end
        else


            behavFile = [dirInput sub '\behav\MSL_openLoop_' sub '_PVT.txt'];

            fid=fopen(behavFile);
            aline = fread(fid, 'char=>char');          % returns a single long string
            fclose(fid);

            aline(aline==uint8(sprintf('\r'))) = [];        % remove cariage return
            aline = tokenize(aline, uint8(sprintf('\n')));  % split on newline

            lines= find(~cellfun(@isempty,regexp(aline, '.*PVT.*')));

            fprintf('%i\n',length(lines)/2)

            lines= reshape(lines,2,length(lines)/2);

            for idx_sess= 1 : size(lines,2)
                sessionLine = strsplit(aline{lines(1,idx_sess)},'-');
                if strcmp(sessionLine{2}(1:2),'PR')
                    session = 'preNight';
                elseif strcmp(sessionLine{2}(1:2),'PO')
                    session = 'postNight';
                end

                counter=1;

                for idx = 1 : 2:lines(2,idx_sess)-lines(1,idx_sess)-1
                    t0=aline{lines(1,idx_sess)+idx};t0=strsplit(t0,' ');
                    t0=str2num(t0{end});
                    t1=aline{lines(1,idx_sess)+idx+1};t1=strsplit(t1,' ');
                    t1=str2num(t1{end});

                    rt = t1-t0;

                    fprintf(file,'%s; %s; %i; %1.4f\n',sub, session,counter,rt);
                    counter=counter+1;

                end
            end
        end

    end
    fclose(file);
end