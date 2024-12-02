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



[listSub,listSubEEG,listSubBehav] = getFullDatasets;
% Donn???es d???crites par la sous-fonction pm_data
data.subjects = listSubBehav;

% indice_effet{1} = 'Sn(X) iRND_SESSIONX_TRAININGxMean^1*bf(1)';


data.All = [1: size(data.subjects,2)]; %% all

data.dir = {[initPath.Exp '\analyses\fmri\Analyses\']};
data.anaffx_date = '031221';

filetot = fopen([initPath.Exp '\analyses\fmri\Analyses\group\task\rfx_' data.anaffx_date '\PrevsPost_betaValues_ActivationAnalysis.csv'],'w');
fprintf(filetot,'%s;%s;%s;%s;%s;%s;\n', 'Sub' , 'ROI_Name','ROI_Coord', 'Session','Condition','Beta');


%------------------------------------------------------------------------

%Pre vs Post Up

input_coord = struct('name',[],'coord',[]);


% NOM Colonnes
%%%%%%%%%%%%%%

indice_effet{1} = 'Sn(1) up_practicepreNightTraining*bf(1)';
indice_effet{2} = 'Sn(1) down_practicepreNightTraining*bf(1)';
indice_effet{3} = 'Sn(1) not_practicepreNightTraining*bf(1)';
indice_effet{4} = 'Sn(3) up_practicepostNightTraining*bf(1)';
indice_effet{5} = 'Sn(3) down_practicepostNightTraining*bf(1)';
indice_effet{6} = 'Sn(3) not_practicepostNightTraining*bf(1)';

input_coord(1).name = 'Putamen_Right'; input_coord(1).coord = [26 -8 -4];
input_coord(2).name = 'M1_Left'; input_coord(2).coord = [-38 -20 52];
input_coord(3).name = 'Hippocampus_Right'; input_coord(3).coord = [32 -38 -6];
input_coord(4).name = 'Caudate_Right'; input_coord(4).coord = [10 0 14];

input_coord(5).name = 'Caudate_R_upVsDown'; input_coord(5).coord = [20 18 12];
input_coord(6).name = 'Caudate_R_upVsNot'; input_coord(6).coord = [18 28 4];
input_coord(7).name = 'Caudate_R_DownVsNot'; input_coord(7).coord = [16 -2 26];
input_coord(8).name = 'Hippocampus_R_upVsNot'; input_coord(8).coord = [36 -36 -4];


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
    
    cd ([char(data.dir),'/',char(data.subjects(data.All(isub))) '/task/ffx/analysisFfx_' char(data.anaffx_date)])
    
    load SPM.mat
    
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
            
            [u s v] = svd(betatmp(:,icol),0);
            % compute the average beta that account for at least 90% of the variance across the 26 voxels
            indiceigen = find(cumsum(diag(s).^2/sum(diag(s).^2)) >=.90); indiceigen = indiceigen(1);
            betatmp2print = mean([u(:,indiceigen)*s(indiceigen,indiceigen)*v(:,indiceigen)'])';
            
            if iU <4
                session = 'pre';
            else
                session = 'post';
            end
            
            if iU == 1 | iU == 4
                condition = 'up';
            elseif iU == 2 | iU == 5
                condition = 'down';
            elseif iU == 3 | iU == 6
                condition = 'not';
            end

            fprintf(filetot,'%s;%s;%s;%s;%s;%1.4f;\n', char(data.subjects(data.All(isub))) , char(input_coord(jcoord).name), num2str(input_coord(jcoord).coord),session,condition,betatmp2print);
            fprintf('Computing %s in ROI %s of subject %s in condition %s\n', session,char(input_coord(jcoord).name),char(data.subjects(data.All(isub))),condition)

        end
    end
end
    
fclose(filetot)
