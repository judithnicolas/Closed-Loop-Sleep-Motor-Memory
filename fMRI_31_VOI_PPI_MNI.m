function fMRI_31_VOI_PPI_MNI(pathIn,listSub,date,voiname,voicoord,contrastname,contrast_nb)
%
% This script computes the VOI and the PPI for PPI analyses in SPM12

% case 1 : VOI
% case 2 : PPI

% Specify contrast of interest, VOIname and VOIcoord
% GA 02/12/2010 Montreal
% MV 21/10/2020 Leuven: add condition factor within session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;
% clc;

global data
spm fmri

warning off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SPECIFY YOUR DATA HERE

data.dir = [pathIn '\analyses\fmri\Analyses\'];

%% Set parameters


% depending on contrast different participants, their corresponding
% individual targets and sessions have to be defined

data.subjects  =   listSub;

% % ------------------------------------------------------------------------
matlabbatch = {};

fprintf(1,'BUILDING JOB: Create VOIs...')
for idx_voi = 1 : length(voiname)
    for isub = 1 : size(data.subjects,2)
        
        cd([data.dir, data.subjects{isub},'/task/ffx/analysisFfx_' date '/raw/']) %from activation based
        load([data.dir, data.subjects{isub},'/task/ffx/analysisFfx_' date '/contrast/contrast_raw.mat']) %from activation based
        
        counter = 1;
        for ises = [1 3]%1 = pre-night training; 2 = post-night training;
            
            
            fprintf(1,'PROCESSING CONTRAST : %s\n',contrastes(contrast_nb(counter)).name)
            matlabbatch{counter}.spm.util.voi.spmmat = {(strcat(data.dir,char(data.subjects{isub}),[ '/task/ffx/analysisFfx_' date '/raw/SPM.mat']))};
            matlabbatch{counter}.spm.util.voi.adjust = 0;
            matlabbatch{counter}.spm.util.voi.session = ises;
            matlabbatch{counter}.spm.util.voi.name = [voiname{idx_voi} '_' contrastname];
            matlabbatch{counter}.spm.util.voi.roi{1}.spm.spmmat = {''};
            matlabbatch{counter}.spm.util.voi.roi{1}.spm.contrast = contrast_nb;
            matlabbatch{counter}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
            matlabbatch{counter}.spm.util.voi.roi{1}.spm.thresh = 0.9; % Threshold for display at the individual level
            matlabbatch{counter}.spm.util.voi.roi{1}.spm.extent = 0;
            matlabbatch{counter}.spm.util.voi.roi{1}.spm.mask = 0;
            
            %when you use a sphere
            matlabbatch{counter}.spm.util.voi.roi{2}.sphere.centre = voicoord{idx_voi};
            matlabbatch{counter}.spm.util.voi.roi{2}.sphere.radius = 6; %cortical = 10; subcortical = 6 ??
            matlabbatch{counter}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
            matlabbatch{counter}.spm.util.voi.expression = 'i1&i2';
            
            %when you use a mask
%                             matlabbatch{1}.spm.util.voi.roi{3}.mask.image = {'D:\KUL\Experiments\Experiment 2 - Neural correlates\Analysis\MRI\Analysis\Mask\Oxford_Mask\harvardoxford-subcortical_prob_Right Hippocampus.nii,1'};
%                             matlabbatch{1}.spm.util.voi.roi{3}.mask.threshold = 0.5;
%                             matlabbatch{1}.spm.util.voi.expression = 'i1&i2&i3'; %when you use a mask

            counter = counter + 1;
        end
        
        eval (['save ' data.dir '/JOBS/jobs_VOI_' voiname{idx_voi} '_' contrastname ' matlabbatch'])
        
        spm_jobman('run',matlabbatch)
        
    end
    
end

fprintf(1,'BUILDING JOB: Create PPIs...')
for idx_voi = 1 : length(voiname)
    
    for isub = 1 : size(data.subjects,2)
        cd([data.dir, data.subjects{isub},'/task/ffx/analysisFfx_' date '/raw/']) %from activation based
        a = dir(strcat('VOI_',voiname{idx_voi},'_',contrastname,'_', num2str(1),'.mat')) %1 refers to training session ??
        if size(a,1)~=0
            for ises = [1 3] % 1 = pre-night training; 2 = post-night training
                matlabbatch = {};
                for icond = 1:3 % 1 = up; 2 = down;  3 = not
                    %add cond
                    matlabbatch{icond}.spm.stats.ppi.spmmat = {(strcat(data.dir,char(data.subjects{isub}), ['/task/ffx/analysisFfx_' date '/raw/SPM.mat']))};
                    matlabbatch{icond}.spm.stats.ppi.type.ppi.voi = {(strcat(data.dir,char(data.subjects{isub}), ['/task/ffx/analysisFfx_' date '/raw/VOI_',voiname{idx_voi},'_'], contrastname,'_', num2str(ises), '.mat'))};
                    
                    % adjust ppi.u based on condition
                    % ppi.u is a [nX3] matrix, with [SPM.Sess.U(i), SPM.Sess.U(i).name{j}, contrast weight]
                    if icond == 1
                        matlabbatch{icond}.spm.stats.ppi.type.ppi.u = [1 1 1];
                    elseif icond == 2
                        matlabbatch{icond}.spm.stats.ppi.type.ppi.u = [2 1 1];
                    elseif icond == 3
                        matlabbatch{icond}.spm.stats.ppi.type.ppi.u = [3 1 1];
                    end
                    
                    %example ppi.u: [1 1 1; 2 1 -1]
                    %first row: first condition (e.g., TMR)
                    %second row: second condition (e.g., Random)
                    %third column: contrast ('1' for the first condition and '-1' for the second condition in the example)
                    %second column: type of analysis ('1' for main effect and '2' for parametric modulation)
                    
                    %Modulation PPI for two conditions -> for each ises:
                    %cond = 1: ppi.u = [1 2 1]
                    %cond = 2: ppi.u = [2 2 1]
                    
                    matlabbatch{icond}.spm.stats.ppi.name = strcat(voiname{idx_voi},'_',contrastname,'_',num2str(ises),'_',num2str(icond));
                    matlabbatch{icond}.spm.stats.ppi.disp = 1; % 1 display, 0 no display
                end %icond
                
                eval (['save ' data.dir '/JOBS/jobs_VOI_' voiname{idx_voi} '_' contrastname ' matlabbatch'])
            
                spm_jobman('run',matlabbatch)
                
            end %ises
            

        end
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%