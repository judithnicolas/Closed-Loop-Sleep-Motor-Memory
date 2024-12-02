function fMRI_11_preprocessing_MNI(path,sub)

%% SPM12 Spatial Preprocessing pipeline - Stress Study - Leuven 12/07/2018 - GA

% Data need to be organized as follows:
% ana\subj_directory\mprage: s*.img
%                   \session1: f*.img 20 blocks
%                   \session2: f*.img 4 blocks
%                   \session3: f*.img 20 blocks

% A good habit is to save the raw data organized as mentioned above in a separate \raw_organized folder.
% We therefore recommend saving the following \raw,\raw_organized and run as many ana as wished in separate \ana folder
% WARNING! ALWAYS RUN NEW ANA ON RAW DATA COPIED FROM RAW_ORGANIZED FOLDER

% INPUT PARAMETERS / SUBJECT INFO

spm_jobman('initcfg');
data.main_folder = {[path '\analyses\fmri\Analyses\']};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% STRESS MSL MR S58826 
%%%% Following participants are not included with dartel: 
%%%%    cortisol responder in control group (n=4) [S17, S24, S46, S47]
%%%%    performance outliers (n=2) [S29, S44]
%%%%    Motion [S21] 
%%%% TOTAL INCLUDED: 73 
% data.subjects   =   getFullDatasets;
data.subjects   = { sub   };
            
for sub = 1:1:length(data.subjects)   

       data.sessions (sub, 1:3) =     { ...
        'prenight/training', ...
        'prenight/test', ...
        'postnight/training',...
    }';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% RE-ORIENT AND SEGMENT T1
% 
% % RE-ORIENT function to automatically reorient a structural image on the template
% % in order further proceed with the segmentation/normalisation.
% % very useful as segmentation/normalisation is rather sensitive on the
% % starting orienation of the image
% 
% % SEGMENT
% % Data classified into different tissue types, typically with the following
% % order: grey matter, white matter, CSF, bone, soft tissue, background.
% 
% % Within each subject's mprage folder, multiple files will be created.
% % Files with prefixes c1-c6 are tissue class images (one for each type of
% % tissue) that are in alignment with the original space (i.e., segmented
% % gray matter, segmented white matter, etc.). Files with the prefixes
% % rc1-rc6 are in a form that can be used with the Dartel toolbox.
% % Forward deformation field ”y_”imagename“_sn.mat” is also written and will
% % be used for the normalize to MNI step (inverse deformation field ”iy_”imagename“_sn.mat”is also
% % written).
% 
o_matlabbatch = [];
for nSub =1:length(data.subjects)
    
        fprintf(['AUTO-ORIENT SUBJECT ', data.subjects{nSub},'\n']);
        anat_folder = fullfile(data.main_folder{1}, data.subjects{nSub},'\task\anat',filesep);
        anat = fMRI_0_get_files(anat_folder,'nii','DBIEX'); % Get nifti files that begin with DBIEX (original)
        %anat = fMRI_0_get_files(anat_folder,'img','s'); % FOR ME TO PLAY
        fMRI_0_reorient(anat);
        fprintf(['DONE AUTO-ORIENT T1 SUBJECT ', data.subjects{nSub},'\n']);
        o_matlabbatch{end+1} = p_segment(anat); % Still retrieving nifti files that begin with DBIEX (but these have not been auto-oriented; not renamed)
    
end

spm_jobman('run',o_matlabbatch);
fprintf(['DONE SEGMENTATION T1 SUBJECT ', data.subjects{nSub},'\n']);

% % %%   REALIGN FUNCTIONAL
% % % Within each subject's primary (functional) data folder, a text file with
% % % the prefix rp will be created with realignment parameters. A nifti file
% % % with the prefix 'mean' will also be created with the mean resliced (functional) image.
% % % No new realigned files (i.e., niftis with prefix 'r' will be written with parameters specified in
% % % p_realign).
 
o_matlabbatch = [];
for nSub=1:length(data.subjects)
    o_matlabbatch{end+1} = p_realign;
    for nRun = 1:size(data.sessions,2)
        realign_folder = fullfile(data.main_folder{1}, data.subjects{nSub},'\task\',data.sessions{nSub,nRun}, filesep);
        realign = fMRI_0_get_files(realign_folder, 'nii','vol');
        realign = strcat(realign, ',1');
        o_matlabbatch{end}{end}.spm.spatial.realign.estwrite.data{nRun} = cellstr(realign);
    end
end
spm_jobman('run',o_matlabbatch);
clear nSub; clear nRun;

%% COREGISTRATION
% Within each subject's primary (anat) data folder, an mstructural file
% will be created (anatomical coregistered). No reslicing of mean and
% functional (estimate only)

o_matlabbatch = [];
for nSub=1:length(data.subjects)
        
        anat_folder = fullfile(data.main_folder{1}, data.subjects{nSub},'\task\anat',filesep);
        anat = fMRI_0_get_files(anat_folder,'nii','DBIEX'); % Get nifti files that begin with DBIEX (original)
        mean_folder = fullfile(data.main_folder{1}, data.subjects{nSub},data.sessions{nSub,1}, filesep);
        mean = fMRI_0_get_files(mean_folder,'nii','mean');
        
        other = [];
        for nRun = 1:size(data.sessions,2)
            realign_folder = fullfile(data.main_folder{1}, data.subjects{nSub},'task\',data.sessions{nSub, nRun}, filesep);
            filesTmp = fMRI_0_get_files(realign_folder, 'nii','vol');
            filesTmp = strcat(filesTmp, ',1');
            other = strvcat(other, filesTmp);
        end
        o_matlabbatch{end+1} = fMRI_0_coregister(anat,mean,other); % anatomical is the reference; mean functional is the source and other is all other functionals are to be converted
     
end
spm_jobman('run',o_matlabbatch);
clear nSub; clear nRun;


%% NORMALIZATION STRUCTURAL AND FUNCTIONAL TO MNI
%% Normalize will be done with DARTEL for the final analysis (i.e., when the full data set is available)

% NORM WRITE ANAT
o_matlabbatch = [];
for nSub=1:length(data.subjects)
        
        anat_folder = fullfile(data.main_folder{1}, data.subjects{nSub},'\task\anat',filesep);
        i_files2normalise = fMRI_0_get_files(anat_folder,'nii','DBIEX'); % Get nifti files that begin with DBIEX (original)
        i_flowfield = fMRI_0_get_files(anat_folder,'nii','y');
        o_matlabbatch{end+1} = fMRI_0_normalise2mni_anat_nodartel(i_flowfield,i_files2normalise);
end
spm_jobman('run',o_matlabbatch);
fprintf(['DONE NORMALIZE ANAT SUBJECT ', data.subjects{nSub},'\n']);
clear nSub;

% Within each subject's anatomical data folder, a new nifti file
% will be created with the prefix 'w' which represents the anatomical scan
% normalized to the MNI template.

% NORM WRITE FUNC
o_matlabbatch = [];
for nSub=1:length(data.subjects)
%     if isempty(dir([data.main_folder{1} '/' data.subjects{nSub} '/' data.sessions{1} '/*.txt']))
        
        i_files2normalise = [];
        for nRun = 1:size(data.sessions,2)
            %         if nRun == 3 && (strcmp(data.subjects{nSub},'S3') || strcmp(data.subjects{nSub},'S4') || strcmp(data.subjects{nSub},'S11')|| strcmp(data.subjects{nSub},'S12') || strcmp(data.subjects{nSub},'S66'))
            %             break; end
            realign_folder = fullfile(data.main_folder{1}, data.subjects{nSub},'\task\',data.sessions{nSub,nRun}, filesep);
            filesTmp = fMRI_0_get_files(realign_folder, 'nii','vol');
            filesTmp = strcat(filesTmp, ',1');
            i_files2normalise = strvcat(i_files2normalise, filesTmp);
        end
        anat_folder = fullfile(data.main_folder{1}, data.subjects{nSub},'\task\anat',filesep);
        i_flowfield = fMRI_0_get_files(anat_folder,'nii','y');
        o_matlabbatch{end+1} = p_normalise2mni_func_nodartel(i_flowfield,i_files2normalise);
%     end
end
spm_jobman('run',o_matlabbatch);
fprintf(['DONE NORMALIZE FUNC SUBJECT ', data.subjects{nSub},'\n']);
clear nSub;

% Within each subject's functional data folder, a new nifti file
% will be created with the prefix 'w' which represents the functional
% scans normalized to the MNI template.

%% SMOOTH FUNCTIONAL IMAGES
%% No smoothing will be applied for the final analysis as NORMALIZE DARTEL includes smoothing
o_matlabbatch = [];
for nSub=1:length(data.subjects)
%     if isempty(dir([data.main_folder{1} '/' data.subjects{nSub} '/' data.sessions{1} '/*.txt']))
        
        i_files2smooth = [];
        for nRun = 1:size(data.sessions,2)
            realign_folder = fullfile(data.main_folder{1}, data.subjects{nSub},'\task\',data.sessions{nSub,nRun}, filesep);
            filesTmp = fMRI_0_get_files(realign_folder, 'nii','wvol');
            filesTmp = strcat(filesTmp, ',1');
            i_files2smooth = strvcat(i_files2smooth, filesTmp);
        end
        o_matlabbatch{end+1} = p_smooth(i_files2smooth);
%     end
end
spm_jobman('run',o_matlabbatch);
fprintf(['DONE SMOOTH FUNC SUBJECT ', data.subjects{nSub},'\n']);
clear nSub;

fprintf('%s','preprocessing done\n')

