function fMRI_32_FFX_PPI_MNI(pathIn,listSub,date,voiname,contrastname)

% GA - Leuven - 24/07/2018 - TO RUN - adjusted by MV to run on TMR-MRI study
% This script  builds up design matrices across a bunch of subjects
% It has to be run in 2 separate steps (action) :
% action 1 = fMRI design and estimate
% action 2 = contrasts
% Your input is twofold :
% 1. Specify your data at the beginning of the script
% 2. Specify your contrasts as a subfunction at the end of the script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global data
spm fmri

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%------------------------------DATA STIMPLAST ------------------------------------%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data.dir = [pathIn '\analyses\fmri\Analyses\'];
data.subjects  =   listSub;

data.ppiname= voiname ;

%%
origdir = pwd;

for idx_voi = 1 : length(voiname)
    matlabbatch = [];
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
    
    for isub = 1  : size(data.subjects,2)
        clear condfile;
        matlabbatch{1}.spm.stats.fmri_spec.sess = {};
        % Select the directory to put SPM.mat
        
        % CHECK ANAFFX FOLDER IN FIRST SESSION
        % IF DOESN'T EXIST CREATE IT OTHERWISE OVERWRITE
        strdir = fullfile(data.dir,char(data.subjects{isub}),['/task/ffx/analysisFfx_' date '/raw/PPI_analysis/' data.ppiname{idx_voi}]);
        if ~exist(char(strdir))
            mkdir(char(strdir))
        else
            fprintf(1,'A DIRECTORY WITH THIS NAME ALREADY EXISTED AND WAS OVERWRITTEN, SORRY \n');
%             rmdir(char(strdir),'s')
            mkdir(char(strdir))
        end
        matlabbatch{1}.spm.stats.fmri_spec.dir = {char(strdir)};
        
        % FMRI DESIGN
        fprintf(1,'BUILDING JOB : FMRI design\n')
        
        % prepare multicondition selection
        % Look for the onsets files in the first directory
        %one onset file in each condition folder, changed by MG:
        %looking for condfile inside loop across conditions
        %??why char(data.sessions{isbu, ses_counter}) instead of char(data.sessions{ses_counter})?
        
        %condsession_order =  [[1:1:size(condfiles,1)]' zeros(size(condfiles,1),1)]; %% when 2 SESSIONS,...
        clear data.sessions
        
        data.sessions(isub,1:2) = {'prenight','postnight'};
        
        data.conditions (isub, 1:3) = { 'up', 'down', 'not'}';
        
        for ises = 1 : size(data.sessions,2) % CONF
            for icond = 1 : size(data.conditions,2)
                % file selection
                [files]=spm_select('List',char([ data.dir '/' char(data.subjects{isub}) '/' char(data.sessions{isub,ises}) '/training/']),'^sw.*\.nii$');
                matlabbatch{1}.spm.stats.fmri_spec.sess(ises).scans =  ...
                    cellstr(strcat(repmat(fullfile(data.dir,char(data.subjects{isub}),char(data.sessions{isub,ises}),'/training/'),size(files,1),1),repmat('/',size(files,1),1),files));
                
                cd(char(strdir))
                cd ../../
                
                if ises ==2
                    sess = 3;
                else
                    sess = ises;
                end
                if isfile(['PPI_' voiname{idx_voi} '_' contrastname '_' num2str(sess) '_' num2str(icond) '.mat'])
                    eval (['load PPI_' voiname{idx_voi} '_' contrastname '_' num2str(sess) '_' num2str(icond)])
                    
                    if icond ==1
                        PPI_ises_cond1 = PPI;
                    elseif icond == 2
                        PPI_ises_cond2 = PPI;
                    elseif icond == 3
                        PPI_ises_cond3 = PPI;
                    end
                end
                    
            end % icond
            
            % Regressors for PPI
            matlabbatch{1}.spm.stats.fmri_spec.sess(ises).regress(1).name = ['Physiological_' data.sessions{isub,ises}];
            matlabbatch{1}.spm.stats.fmri_spec.sess(ises).regress(2).name = ['Psychological_' data.sessions{isub,ises} '_' data.conditions{isub,1}]; %Psych cond1
            matlabbatch{1}.spm.stats.fmri_spec.sess(ises).regress(3).name = ['Psychological_' data.sessions{isub,ises} '_' data.conditions{isub,2}]; %Psych cond2
            matlabbatch{1}.spm.stats.fmri_spec.sess(ises).regress(4).name = ['Psychological_' data.sessions{isub,ises} '_' data.conditions{isub,3}]; %Psych cond2
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(ises).regress(5).name = ['PPI_' data.sessions{isub,ises} '_' data.conditions{isub,1}];% PPI cond1
            matlabbatch{1}.spm.stats.fmri_spec.sess(ises).regress(6).name = ['PPI_' data.sessions{isub,ises} '_' data.conditions{isub,2}];% PPI cond2
            matlabbatch{1}.spm.stats.fmri_spec.sess(ises).regress(7).name = ['PPI_' data.sessions{isub,ises} '_' data.conditions{isub,3}];% PPI cond2
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(ises).regress(1).val = PPI.Y; % PPI_isess_icond1_Y
            matlabbatch{1}.spm.stats.fmri_spec.sess(ises).regress(2).val = PPI_ises_cond1.P; % PPI cond 1.P
            matlabbatch{1}.spm.stats.fmri_spec.sess(ises).regress(3).val = PPI_ises_cond2.P; % PPI cond 2.P
            matlabbatch{1}.spm.stats.fmri_spec.sess(ises).regress(4).val = PPI_ises_cond3.P; % PPI cond 2.P
            matlabbatch{1}.spm.stats.fmri_spec.sess(ises).regress(5).val = PPI_ises_cond1.ppi; % PPI cond 1.ppi
            matlabbatch{1}.spm.stats.fmri_spec.sess(ises).regress(6).val = PPI_ises_cond2.ppi; % PPI cond 2.ppi
            matlabbatch{1}.spm.stats.fmri_spec.sess(ises).regress(7).val = PPI_ises_cond3.ppi; % PPI cond 2.ppi
            
            % multiregressor selection
            rpfile=spm_select('List',char(fullfile(data.dir,char(data.subjects{isub}),char(data.sessions{isub,ises}),'training/')),strcat('^rp.*\.txt$'));
            matlabbatch{1}.spm.stats.fmri_spec.sess(ises).multi_reg = ...
                cellstr(fullfile(data.dir,char(data.subjects{isub}),char(data.sessions{isub,ises}),'training/',rpfile));

            % HPF
            matlabbatch{1}.spm.stats.fmri_spec.sess(ises).hpf = 128;
        end % ises
        
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        
        
        % FMRI ESTIMATE
        fprintf(1,'BUILDING JOB : FMRI estimate\n')
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname = 'fMRI model specification: SPM.mat File';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        
        
        cd(char(strdir))
        eval (['save jobs_TMR-MRI_FFX_PPI_SPM12_' data.subjects{isub} '_' data.ppiname{idx_voi} ' matlabbatch'])
        spm_jobman('run',matlabbatch)
        contrastes = get_contrast_for_PPI_ffx(strdir);

    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'BUILDING JOB : FMRI contrasts\n')
for idx_voi = 1 : length(voiname)
    
    for isub = 1  :size(data.subjects,2)
        
        load (fullfile(data.dir,char(data.subjects{isub}),['/task/ffx/analysisFfx_' date '/raw/PPI_analysis/contrastes.mat' ]))
        
        fprintf(1,'PROCESSING SUBJECT %i : %s\n',isub,data.subjects{isub})
        
        strdir = fullfile(data.dir,char(data.subjects{isub}),['/task/ffx/analysisFfx_' date '/raw/PPI_analysis/'],voiname{idx_voi});
        cd(strdir)
        eval (['load jobs_TMR-MRI_FFX_PPI_SPM12_' data.subjects{isub} '_' voiname{idx_voi}  '.mat']) % Loading the ss's job to append the contrasts
        
        matlabbatch{3} = {};
        matlabbatch{3}.spm.stats.con.spmmat = cellstr(fullfile(strdir,'SPM.mat'));
        
        for icon = 1:size(contrastes,2)
            matlabbatch{3}.spm.stats.con.consess{icon}.tcon.name = contrastes(icon).name;
            matlabbatch{3}.spm.stats.con.consess{icon}.tcon.convec = contrastes(icon).C;
            matlabbatch{3}.spm.stats.con.consess{icon}.tcon.sessrep = 'none';
        end % icon
        
        matlabbatch{3}.spm.stats.con.delete = 0;
        
        eval (['save jobs_TMR-MRI_FFX_PPI_SPM12_' data.subjects{isub} '_' data.ppiname{idx_voi} '.mat ' ' matlabbatch'])
        
        % Pruning the empty options
        for ii = 1:2
            matlabbatch(1) = []; % no curly brackets to remove array
        end
        
        spm_jobman('run',matlabbatch)
        
        
        clear matlabbatch
    end % isub dans contrastes
    
end
cd(char(data.dir))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
