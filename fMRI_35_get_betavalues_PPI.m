% Collect betas for target voxels over subjects
% Look for local maxima within 3mm from the target voxel
% From the target voxel, further look around it there is a maximum in the direct neighbours
% Compute mean and std betas over subjects
% Display results for each target area

% TO BE ENTERED
% coordaentrer : specifying the name and coordinates (mm) of the voxel
% indices : les betas d'interet (colonnes de X)
% data  : the directories of the various subjects (same as in pm_sptr*, pm_ffx*, or pm_rfx*)
%
% pm 8/6/4
% ---------------------------------------------------------------------------------------------
% %  Transforme les coordonn�es voxels and mm
% XYZmm = SPM.xVol.M(1:3,:)*[SPM.xVol.XYZ; ones(1,size(SPM.xVol.XYZ,2))]
% %  Transforme les coordonn�es mm en voxels
% inv(SPM.xVol.M)*[kk;ones(1,size(kk,2))]
% ---------------------------------------------------------------------------------------------

% ENTREES
% Entrer les coordonn�es � rechercher (plusieurss coordonn�es NE peuvent PAS etre trait�es simultan�ment)
%Changer le nom de la .ma cree premier item apres save de la ligne 243.


[TMRup,TMRdown,offlineGainUp,offlineGainDown,offlineGainNot]= getBehaviour();

[listSub,listSubEEG,listSubBehav] = getFullDatasets;
data.subjects = listSubBehav;

data.All = [1: size(data.subjects,2)]; %% all

data.dir = {[initPath.Exp '\analyses\fmri\Analyses\']};
data.anaffx_date = '031221';


%------------------------------------------------------------------------


%Pre vs Post Up

input_coord = struct('name',[],'coord',[]);


% NOM Colonnes
%%%%%%%%%%%%%%

indice_effet{1} = 'Sn(1) PPI_prenight_up';
indice_effet{2} = 'Sn(1) PPI_prenight_down';
indice_effet{3} = 'Sn(1) PPI_prenight_not';
indice_effet{4} = 'Sn(2) PPI_postnight_up';
indice_effet{5} = 'Sn(2) PPI_postnight_down';
indice_effet{6} = 'Sn(2) PPI_postnight_not';


%% Right Caudate------------------------------------------------------------------------


seed = 'right_caudate';

filetot = fopen([initPath.Exp '\analyses\fmri\Analyses\group\task\rfx_' data.anaffx_date '\PPI_' seed '.csv'],'w');
fprintf(filetot,'%s;%s;%s;%s;%s;%s;%s;%s;%s;\n', 'Sub' , 'ROI_Name','ROI_Coord','Condition', 'Session','Contrast','Beta');



%Pre vs Post Up

input_coord = struct('name',[],'coord',[]);


input_coord(1).name = 'hippocampus_left';   input_coord(1).coord = [-32 -20 -14]; 


input_contrast{1} = 'overnightChange_Not';     

%-------------------------------------------------------------------------------------
origdir = pwd;
BETA = [];% BETA est une matrice nbeta x ncoord x nsujet
seedcoordmm = [];
for dd = 1:size(input_coord,2)
    seedcoordmm = [seedcoordmm input_coord(dd).coord'];
end
shift = [repmat([-1 0 1]',9,1) ...
    repmat(kron([-1 0 1]',ones(3,1)),3,1) ...
    repmat(kron([-1 0 1]',ones(9,1)),1,1)];

%------------------------------------------------------------------------------------
% GET BETA
%-------------------------------------------------------------------------------------
% get betas over subjects for every area
% FORMAT spm_progress_bar('Init',height,xlabel,ylabel,flgs)
% Initialises the bar in the 'Interactive' window.
% If flgs contains a 't', then use tex interpreter for labels.
%
% FORMAT spm_progress_bar('Set',value)
% Sets the height of the bar itself.
%
% FORMAT spm_progress_bar('Clear')
% Clears the 'Interactive' window.



%%% ALL
%%%%%
%%%%%
BETA = [];
RES = [];
%spm_progress_bar('Init',10,1,1,'t')
coord_used = [];
for isub = 1 : size(data.All,2)
    betacoord = [];
    
    cd ([char(data.dir),'/',char(data.subjects(data.All(isub))) '/task/ffx/analysisFfx_' char(data.anaffx_date) '/raw/PPI_analysis/' seed '/'])
    
    load SPM.mat
    counter = 1;

    coord = [];
    seedcoord = inv(SPM.xVol.M)*[seedcoordmm;ones(1,size(seedcoordmm,2))];
    % Search for the local maximum
    for jcoord = 1 : size(seedcoord,2) % pour chaque ROI
        
        betatmp = [];
        betatmp2 = [];
        targetcoord = [];
        
        [xyz,i,d] = spm_XYZreg('NearestXYZ',seedcoord(:,jcoord),SPM.xVol.XYZ); % coordonnees en voxels
        xyzmm = round(SPM.xVol.M(1:3,:)*[xyz; ones(1,size(xyz,2))]);  % nouvelles coordonn�es en mm
        if sqrt(sum((seedcoordmm(1:3,jcoord)-xyzmm).^2)) > 3 % si le target voxel est � plus de 3 mm du seed voxel
            targetcoord = seedcoord(1:3,jcoord);
        else
            targetcoord =  xyz;
        end
        coord = [coord; targetcoord];
        % Search for the max (mean) beta value around the target voxel for each ROI
        for ivox = 1:27 % voxel autour du voxel cible
            tmp = spm_get_data(SPM.Vbeta, targetcoord+shift(ivox,:)');
            betatmp = [betatmp;tmp'];
        end
        % Check that all effects of interest exist in this particular
        % subject and label any absent event type by zeros
        
        % MG: skip first characters SN(?) and SESSIONSX use cmpstr, e.g.
        % strcmp(SPM.xX.name{ii}(7:10), 'cMSL') &&strcmp(SPM.xX.name{ii}(21:end), 'TRAINING*bf(1)')
        for iU = 1 : size(indice_effet,2)
            Xcolnum = [];TEST = [];
            for icol = 1:size(SPM.xX.X,2)
                %                 if length(SPM.xX.name{icol}) > 21
                test = findstr(char(indice_effet{iU}),char(SPM.xX.name(icol)));
                %                     test = strcmp(indice_effet{iU}(6:10),SPM.xX.name{icol}(6:10))&&strcmp(indice_effet{iU}(21:end),SPM.xX.name{icol}(21:end)); %compare cond(cMSL,etc) and endpar (Training vs TrainingxMEAN
                if test
                    betatmp2 = [betatmp2 betatmp(:,icol)];
                    TEST = [TEST 1];
                    break
                else
                    TEST = [TEST 0];
                end
                %                 else
                %                     TEST = [TEST 0];
                %                 end
            end
            if ~any(TEST)
                betatmp2 = [betatmp2 zeros(size(betatmp,1),1)];
                betatmp2(isnan(betatmp(:,end-1)),end) = NaN; % ne remplacer que les valeurs num�riques sinon le reshape plus bas ne marche plus
            end
            
            [u s v] = svd(betatmp(~isnan(betatmp(:,icol)),icol),0);
            % compute the average beta that account for at least 90% of the variance across the 26 voxels
            indiceigen = find(cumsum(diag(s).^2/sum(diag(s).^2)) >=.90); indiceigen = indiceigen(1);
            betatmp2print = mean([u(:,indiceigen)*s(indiceigen,indiceigen)*v(:,indiceigen)'])';
            
            counter = counter+1;
            if iU ==1 | iU == 2 | iU == 3
                session = 'pre';
            else
                session = 'post';
            end
            
            if iU ==1 | iU ==4
                cond = 'up';
            elseif iU ==2 | iU ==5
                cond = 'down';
            elseif iU ==3 | iU ==6
                cond = 'not';
            end
            
            fprintf(filetot,'%s;%s;%s;%s;%s;%s;%1.4f\n', char(data.subjects(data.All(isub))) , ...
                char(input_coord(jcoord).name), num2str(input_coord(jcoord).coord),...
                cond,session,input_contrast{jcoord},betatmp2print);
            fprintf('Computing %s in ROI %s of subject %s in condition %s  beta : %1.4f\n', ...
                session,num2str(input_coord(jcoord).coord),char(data.subjects(data.All(isub))),indice_effet{iU},betatmp2print)

        end
    end
end
    

fclose(filetot);


%% Right Hippocampus------------------------------------------------------------------------
seed = 'right_hippocampus';

filetot = fopen([initPath.Exp '\analyses\fmri\Analyses\group\task\rfx_' data.anaffx_date '\PPI_' seed '.csv'],'w');
fprintf(filetot,'%s;%s;%s;%s;%s;%s;\n', 'Sub' , 'ROI_Name','ROI_Coord','Condition', 'Session','Beta');


input_coord(1).name = 'pmc_right';   input_coord(1).coord = [46 -8 54]; 
input_coord(2).name = 'aspl_right';   input_coord(2).coord = [60 -24 48]; 
input_coord(3).name = 'M1_right_downVsNot';   input_coord(3).coord = [54 -22 46]; % same for up vs not and down vs not
input_coord(4).name = 'aspl_right_shifted';   input_coord(4).coord = [58 -22 46];
input_coord(4).name = 'right_M1';   input_coord(4).coord = [52 -24 40];

 
%-------------------------------------------------------------------------------------
origdir = pwd;
BETA = [];% BETA est une matrice nbeta x ncoord x nsujet
seedcoordmm = [];
for dd = 1:size(input_coord,2)
    seedcoordmm = [seedcoordmm input_coord(dd).coord'];
end
shift = [repmat([-1 0 1]',9,1) ...
    repmat(kron([-1 0 1]',ones(3,1)),3,1) ...
    repmat(kron([-1 0 1]',ones(9,1)),1,1)];

%------------------------------------------------------------------------------------
% GET BETA
%-------------------------------------------------------------------------------------
% get betas over subjects for every area
% FORMAT spm_progress_bar('Init',height,xlabel,ylabel,flgs)
% Initialises the bar in the 'Interactive' window.
% If flgs contains a 't', then use tex interpreter for labels.
%
% FORMAT spm_progress_bar('Set',value)
% Sets the height of the bar itself.
%
% FORMAT spm_progress_bar('Clear')
% Clears the 'Interactive' window.



%%% ALL
%%%%%
%%%%%
BETA = [];
RES = [];
%spm_progress_bar('Init',10,1,1,'t')
coord_used = [];
for isub = 1 : size(data.All,2)
    betacoord = [];
    
    cd ([char(data.dir),'/',char(data.subjects(data.All(isub))) '/task/ffx/analysisFfx_' char(data.anaffx_date) '/raw/PPI_analysis/' seed '/'])
    
    load SPM.mat
    counter = 1;

    coord = [];
    seedcoord = inv(SPM.xVol.M)*[seedcoordmm;ones(1,size(seedcoordmm,2))];
    % Search for the local maximum
    for jcoord = 1 : size(seedcoord,2) % pour chaque ROI
        
        betatmp = [];
        betatmp2 = [];
        targetcoord = [];
        
        [xyz,i,d] = spm_XYZreg('NearestXYZ',seedcoord(:,jcoord),SPM.xVol.XYZ); % coordonnees en voxels
        xyzmm = round(SPM.xVol.M(1:3,:)*[xyz; ones(1,size(xyz,2))]);  % nouvelles coordonn�es en mm
        if sqrt(sum((seedcoordmm(1:3,jcoord)-xyzmm).^2)) > 3 % si le target voxel est � plus de 3 mm du seed voxel
            targetcoord = seedcoord(1:3,jcoord);
        else
            targetcoord =  xyz;
        end
        coord = [coord; targetcoord];
        % Search for the max (mean) beta value around the target voxel for each ROI
        for ivox = 1:27 % voxel autour du voxel cible
            tmp = spm_get_data(SPM.Vbeta, targetcoord+shift(ivox,:)');
            betatmp = [betatmp;tmp'];
        end
        % Check that all effects of interest exist in this particular
        % subject and label any absent event type by zeros
        
        % MG: skip first characters SN(?) and SESSIONSX use cmpstr, e.g.
        % strcmp(SPM.xX.name{ii}(7:10), 'cMSL') &&strcmp(SPM.xX.name{ii}(21:end), 'TRAINING*bf(1)')
        for iU = 1 : size(indice_effet,2)
            Xcolnum = [];TEST = [];
            for icol = 1:size(SPM.xX.X,2)
                %                 if length(SPM.xX.name{icol}) > 21
                test = findstr(char(indice_effet{iU}),char(SPM.xX.name(icol)));
                %                     test = strcmp(indice_effet{iU}(6:10),SPM.xX.name{icol}(6:10))&&strcmp(indice_effet{iU}(21:end),SPM.xX.name{icol}(21:end)); %compare cond(cMSL,etc) and endpar (Training vs TrainingxMEAN
                if test
                    betatmp2 = [betatmp2 betatmp(:,icol)];
                    TEST = [TEST 1];
                    break
                else
                    TEST = [TEST 0];
                end
                %                 else
                %                     TEST = [TEST 0];
                %                 end
            end
            if ~any(TEST)
                betatmp2 = [betatmp2 zeros(size(betatmp,1),1)];
                betatmp2(isnan(betatmp(:,end-1)),end) = NaN; % ne remplacer que les valeurs num�riques sinon le reshape plus bas ne marche plus
            end
            
            [u s v] = svd(betatmp(~isnan(betatmp(:,icol)),icol),0);
            % compute the average beta that account for at least 90% of the variance across the 26 voxels
            indiceigen = find(cumsum(diag(s).^2/sum(diag(s).^2)) >=.90); indiceigen = indiceigen(1);
            betatmp2print = mean([u(:,indiceigen)*s(indiceigen,indiceigen)*v(:,indiceigen)'])';
            
            counter = counter+1;
            if iU ==1 | iU == 2 | iU == 3
                session = 'pre';
            else
                session = 'post';
            end
            
            if iU ==1 | iU ==4
                cond = 'up';
            elseif iU ==2 | iU ==5
                cond = 'down';
            elseif iU ==3 | iU ==6
                cond = 'not';
            end
            
            fprintf(filetot,'%s;%s;%s;%s;%s;%1.4f\n', char(data.subjects(data.All(isub))) , ...
                char(input_coord(jcoord).name), num2str(input_coord(jcoord).coord),...
                cond,session,betatmp2print);
            fprintf('Computing %s in ROI %s of subject %s in condition %s  beta : %1.4f\n', ...
                session,num2str(input_coord(jcoord).coord),char(data.subjects(data.All(isub))),indice_effet{iU},betatmp2print)

        end
    end
end
    

fclose(filetot);




%% Right Putamen------------------------------------------------------------------------

seed = 'right_putamen';

filetot = fopen([initPath.Exp '\analyses\fmri\Analyses\group\task\rfx_' data.anaffx_date '\PPI_' seed '.csv'],'w');
fprintf(filetot,'%s;%s;%s;%s;%s;%s;\n', 'Sub' , 'ROI_Name','ROI_Coord','Condition', 'Session','Beta');


%Pre vs Post Up

input_coord = struct('name',[],'coord',[]);

input_coord(1).name = 'precentral_right';   input_coord(1).coord = [42 -4 56]; 
input_coord(2).name = 'precentral_left';   input_coord(2).coord = [-32 -20 50]; 
input_coord(3).name = 'M1_right';   input_coord(3).coord = [28 -8 44];

 
%-------------------------------------------------------------------------------------
origdir = pwd;
BETA = [];% BETA est une matrice nbeta x ncoord x nsujet
seedcoordmm = [];
for dd = 1:size(input_coord,2)
    seedcoordmm = [seedcoordmm input_coord(dd).coord'];
end
shift = [repmat([-1 0 1]',9,1) ...
    repmat(kron([-1 0 1]',ones(3,1)),3,1) ...
    repmat(kron([-1 0 1]',ones(9,1)),1,1)];

%------------------------------------------------------------------------------------
% GET BETA
%-------------------------------------------------------------------------------------
% get betas over subjects for every area
% FORMAT spm_progress_bar('Init',height,xlabel,ylabel,flgs)
% Initialises the bar in the 'Interactive' window.
% If flgs contains a 't', then use tex interpreter for labels.
%
% FORMAT spm_progress_bar('Set',value)
% Sets the height of the bar itself.
%
% FORMAT spm_progress_bar('Clear')
% Clears the 'Interactive' window.



%%% ALL
%%%%%
%%%%%
BETA = [];
RES = [];
%spm_progress_bar('Init',10,1,1,'t')
coord_used = [];
for isub = 1 : size(data.All,2)
    betacoord = [];
    
    cd ([char(data.dir),'/',char(data.subjects(data.All(isub))) '/task/ffx/analysisFfx_' char(data.anaffx_date) '/raw/PPI_analysis/' seed '/'])
    
    load SPM.mat
    counter = 1;

    coord = [];
    seedcoord = inv(SPM.xVol.M)*[seedcoordmm;ones(1,size(seedcoordmm,2))];
    % Search for the local maximum
    for jcoord = 1 : size(seedcoord,2) % pour chaque ROI
        
        betatmp = [];
        betatmp2 = [];
        targetcoord = [];
        
        [xyz,i,d] = spm_XYZreg('NearestXYZ',seedcoord(:,jcoord),SPM.xVol.XYZ); % coordonnees en voxels
        xyzmm = round(SPM.xVol.M(1:3,:)*[xyz; ones(1,size(xyz,2))]);  % nouvelles coordonn�es en mm
        if sqrt(sum((seedcoordmm(1:3,jcoord)-xyzmm).^2)) > 3 % si le target voxel est � plus de 3 mm du seed voxel
            targetcoord = seedcoord(1:3,jcoord);
        else
            targetcoord =  xyz;
        end
        coord = [coord; targetcoord];
        % Search for the max (mean) beta value around the target voxel for each ROI
        for ivox = 1:27 % voxel autour du voxel cible
            tmp = spm_get_data(SPM.Vbeta, targetcoord+shift(ivox,:)');
            betatmp = [betatmp;tmp'];
        end
        % Check that all effects of interest exist in this particular
        % subject and label any absent event type by zeros
        
        % MG: skip first characters SN(?) and SESSIONSX use cmpstr, e.g.
        % strcmp(SPM.xX.name{ii}(7:10), 'cMSL') &&strcmp(SPM.xX.name{ii}(21:end), 'TRAINING*bf(1)')
        for iU = 1 : size(indice_effet,2)
            Xcolnum = [];TEST = [];
            for icol = 1:size(SPM.xX.X,2)
                %                 if length(SPM.xX.name{icol}) > 21
                test = findstr(char(indice_effet{iU}),char(SPM.xX.name(icol)));
                %                     test = strcmp(indice_effet{iU}(6:10),SPM.xX.name{icol}(6:10))&&strcmp(indice_effet{iU}(21:end),SPM.xX.name{icol}(21:end)); %compare cond(cMSL,etc) and endpar (Training vs TrainingxMEAN
                if test
                    betatmp2 = [betatmp2 betatmp(:,icol)];
                    TEST = [TEST 1];
                    break
                else
                    TEST = [TEST 0];
                end
                %                 else
                %                     TEST = [TEST 0];
                %                 end
            end
            if ~any(TEST)
                betatmp2 = [betatmp2 zeros(size(betatmp,1),1)];
                betatmp2(isnan(betatmp(:,end-1)),end) = NaN; % ne remplacer que les valeurs num�riques sinon le reshape plus bas ne marche plus
            end
            
            [u s v] = svd(betatmp(~isnan(betatmp(:,icol)),icol),0);
            % compute the average beta that account for at least 90% of the variance across the 26 voxels
            indiceigen = find(cumsum(diag(s).^2/sum(diag(s).^2)) >=.90); indiceigen = indiceigen(1);
            betatmp2print = mean([u(:,indiceigen)*s(indiceigen,indiceigen)*v(:,indiceigen)'])';
            
            counter = counter+1;
            if iU ==1 | iU == 2 | iU == 3
                session = 'pre';
            else
                session = 'post';
            end
            
            if iU ==1 | iU ==4
                cond = 'up';
            elseif iU ==2 | iU ==5
                cond = 'down';
            elseif iU ==3 | iU ==6
                cond = 'not';
            end
            
            fprintf(filetot,'%s;%s;%s;%s;%s;%1.4f\n', char(data.subjects(data.All(isub))) , ...
                char(input_coord(jcoord).name), num2str(input_coord(jcoord).coord),...
                cond,session,betatmp2print);
            fprintf('Computing %s in ROI %s of subject %s in condition %s  beta : %1.4f\n', ...
                session,num2str(input_coord(jcoord).coord),char(data.subjects(data.All(isub))),indice_effet{iU},betatmp2print)

        end
    end
end
    

fclose(filetot);

