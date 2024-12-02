function fMRI_0_smoothAllContrastImages(pathIn,listSubBehav,date,mod,extension)

% This script smooth all con images created at the fisrt level in each
% subject, create a mean structural image and mean mask over the
% population, process the factorial design specification (one sample t
% test) and estimation and estimate Contrats.

% action 1 : smooth con images DO NOT SMOOTH IF BEEN SMOOTHED IN RFX ONE
% SAMPLE
% action 2 : Mean Struct, MeanMask, Factorial design specification and
% estimation, Contrast estimation

% You input is twofold :
% 1. Specify your data at the beginning of the script
% 2. In the origdir, create a structure containing all possible contrasts,
% example: ConOfInterest.mat containing structures Training.con (commun to all subjects)
% and Session.con (contrasts of interest by group)

% ga 17/09/2009 Montreal
% fg 26/03/2012
% mg 23/11/2018 Leuven
% jn 031221 SLC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spm fmri

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%------------------------------DATA STIMPLAST ------------------------------------%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod == 0 
    name = 'raw';
else
    name = 'mod';
end


data.dir = [pathIn '\analyses\fmri\Analyses\'];
data.subjects   =   listSubBehav;
data.anatDir = '\anat\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter = 0;
matlabbatch = {};


for idx_sub = 1 : size(data.subjects,2)
    ConFile = dir([data.dir char(data.subjects{idx_sub}) '\task\ffx\analysisFfx_' date '\' name extension '\con_*.nii']);
    for idx_contrastFile = 1:size(ConFile,1)
        counter = counter + 1;
        matlabbatch{1}.spm.spatial.smooth.data{counter,:} = [data.dir char(data.subjects{idx_sub}) '\task\ffx\analysisFfx_' date '\' name extension '\' ConFile(idx_contrastFile).name];
    end
end

matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';

eval (['save ' data.dir '\JOBS\jobs_SmoothCon_SPM12_CLTMR_MNI' ' matlabbatch'])

fprintf(1,'SMOOTHING CON IMAGES...')
spm_jobman('run',matlabbatch)