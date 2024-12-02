function fMRI_13_model_specification_estimation(pathIn,sub,date)

% JN - Salt Lake City - 24/11/2021
% This script  builds up design matrices across a bunch of subjects
% It has to be run in 2 separate steps (action) :
% action 1 = fMRI design and estimate
% action 2 = contrasts
% Your input is twofold :
% 1. Specify your data at the beginning of the script
% 2. Specify your contrasts as a subfunction at the end of the script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global data
spm fmri


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%------------------------------DATA STIMPLAST ------------------------------------%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data.dir = {[pathIn '\analyses\fmri\Analyses\']};

data.subjects   = { sub   };

for sub = 1:1:length(data.subjects)
    
    data.sessions (sub, 1:3) =     { ...
        'prenight/training', ...
        'prenight/test', ...
        'postnight/training',...
        }';
    
end

data.outputName = ['ffx_' char(date) '.mat'];  % Batch will be saved in the subject directory
data.anatDir = '\task\anat\';
data.anaffx_date = ['analysisFfx_' char(date) '.mat'];
mkdir([char(fullfile(data.dir, '\JOBS\')), data.outputName(1:end-4)]);
data.jobs= [ char(data.dir) '\JOBS\' char(data.outputName(1:end-4)) '/'] ;

cd (char(data.dir))
tic
matlabbatch = [];
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
for isub = 1:size(data.subjects,2)%2, not 1 MG
    clear condfile;
    matlabbatch{1}.spm.stats.fmri_spec.sess = {};
    
    % Select the directory to put SPM.mat
    toc
    % CHECK ANAFFX FOLDER IN FIRST SESSION
    % IF DOESN'T EXIST CREATE IT OTHERWISE OVERWRITE
    strdir = fullfile(char(data.dir),char(data.subjects{isub}),'\task\ffx\',char(data.anaffx_date), '/');
    if ~exist(char(strdir))
        mkdir(char(strdir))
    else
        fprintf(1,'A DIRECTORY WITH THIS NAME ALREADY EXISTED AND WAS OVERWRITTEN, SORRY \n');
        rmdir(char(strdir),'s')
        mkdir(char(strdir))
    end
    matlabbatch{1}.spm.stats.fmri_spec.dir = {char(strdir)};
    
    % FMRI DESIGN
    fprintf(1,'BUILDING JOB : FMRI design\n')
    toc
    for idx_sess = 1:size(data.sessions,2)
        if idx_sess == 1
            nbrBlocks = [1 : 21];
            dataFileName = 'preNightTraining';
        elseif idx_sess == 2
            nbrBlocks = [22 : 24];
            dataFileName = 'preNightTest';
        elseif idx_sess == 3
            nbrBlocks = [25 : 45];
            dataFileName = 'postNightTraining';
        end
        % file selection
        [files]=spm_select('List',[char(data.dir ) char(data.subjects{isub}) '\' char(data.sessions(sub, idx_sess)) '\sw*.nii']);
        matlabbatch{1}.spm.stats.fmri_spec.sess(idx_sess).scans = cellstr(strcat(repmat(fullfile(data.dir,char(data.subjects{isub}),char(data.sessions(sub, idx_sess))),size(files,1),1),repmat('/',size(files,1),1),files));
        
        % multicondition selection
        matlabbatch{1}.spm.stats.fmri_spec.sess(idx_sess).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(idx_sess).multi = cellstr([char(data.dir) char( data.subjects{isub}) '\task\onsets\' char( data.subjects{isub}) '_RegressorBehav_' dataFileName '.mat']);
        
        % multiregressor selection
        matlabbatch{1}.spm.stats.fmri_spec.sess(idx_sess).regress = struct('name', {}, 'val', {});
        [regFile]=spm_select('List',[char(data.dir) char( data.subjects{isub}) '/' char(data.sessions(sub, idx_sess)) '/*rp*.txt' ]);
        matlabbatch{1}.spm.stats.fmri_spec.sess(idx_sess).multi_reg = cellstr([char(data.dir) char( data.subjects{isub}) '/' char(data.sessions(sub, idx_sess)) '/' regFile ]);
        
        % HPF
        matlabbatch{1}.spm.stats.fmri_spec.sess(idx_sess).hpf = 128;
    end
    
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FMRI ESTIMATE
    fprintf(1,'BUILDING JOB : FMRI estimate\n')
    toc
    str_ssname = char(data.sessions(sub, 1));
    % PATH TO SPM.MAT
    mat_dir = char(fullfile(data.dir,char(data.subjects{isub}),'\task\ffx\', data.anaffx_date,'\SPM.mat'));
    
    matlabbatch{2}.spm.stats.fmri_est.spmmat = {mat_dir};
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    eval (['save ' data.jobs 'jobs_modelEstimation_CLTMR_MNI' char(fullfile(data.subjects)) '.mat ' ' matlabbatch'])

    % SAVE MATLABBATCH
    %             RUN SPM JOBMAN
    spm_jobman('run',matlabbatch)
    
    mkdir([char(data.dir) char(data.subjects) '\task\modelEstimate\']);

    copyfile (mat_dir, [char(data.dir) char(data.subjects) '\task\modelEstimate\'])

end
