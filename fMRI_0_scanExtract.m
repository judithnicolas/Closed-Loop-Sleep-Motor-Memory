function scanExtract(path,sub)
    
    fid=fopen([ path 'data\' sub '\behav\MSL_openLoop_' sub '.txt']);
    aline = fread(fid, inf, 'char=>char');          % returns a single long string
    fclose(fid);
    
    datafilepointer1 = fopen([ path 'data\' sub '\fmri\cutScans_' sub '.txt'],'wt');
    datafilepointer2 = fopen([ path '\analyses\fmri\RawOrganized\' sub '\cutScans_' sub '.txt'],'wt');
    
    aline(aline==uint8(sprintf('\r'))) = [];        % remove cariage return
    aline = tokenize(aline, uint8(sprintf('\n')));  % split on newline
    
    startScan= find(~cellfun(@isempty,regexp(aline, '.*StartScan.*')));
    
    if length(startScan)>5
        warning('Check text file for multiple sessions')
    end
    
    for idx_MSL = 1 : length(startScan)
        
        task = strsplit(aline{startScan(idx_MSL)-1},'\t');
        task = strsplit(task{end},' ');
        
        if strcmp(task{1}(1:2),'RC')
            
        else


            if strcmp(task{1},'SL-TRAININGPRE_NIGHT')
                stopScan = find(~cellfun(@isempty,regexp(aline, '.*SL-TRAININGPRE_NIGHT STOP.*')));
            elseif strcmp(task{1},'SL-TEST-PRE_NIGHT')
                stopScan = find(~cellfun(@isempty,regexp(aline, '.*SL-TEST-PRE_NIGHT STOP.*')));
            elseif strcmp(task{1},'SL-TRAINING-POST_NIGHT')
                stopScan = find(~cellfun(@isempty,regexp(aline, '.*SL-TRAINING-POST_NIGHT STOP.*')));
            end

            startPractice= find(~cellfun(@isempty,regexp(aline, '.*Practice .*')));
            startPractice= startPractice(startPractice(:)>startScan(idx_MSL) & startPractice(:)<stopScan  );


            tStart= strsplit(aline{startScan(idx_MSL  )},':');
            tStart = str2num(tStart{1});

            endPractice= find(~cellfun(@isempty,regexp(aline, '.*rest.*')));
            endPractice = endPractice(endPractice(:)>startPractice(1) & endPractice(:)<stopScan);
            tEnd= strsplit(aline{endPractice(end)},' ');
            tEnd = str2num(tEnd{end})+10;

            nbScans = floor((tEnd - tStart)/2)+1;

            fprintf(datafilepointer1, '%i\n', nbScans);
            fprintf(datafilepointer2, '%i\n', nbScans);

            fprintf('%i\n', nbScans)
        end
    end
    fclose(datafilepointer1);
    fclose(datafilepointer2);
    fprintf('%s','scans extracted\n')
end
