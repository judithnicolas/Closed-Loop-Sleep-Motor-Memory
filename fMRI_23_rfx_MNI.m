function fMRI_23_rfx_MNI(pathIn,listSubBehav,date,mod,extension)

% This script smooth all con images created at the fisrt level in each
% subject, create a mean structural image and mean mask over the
% population, process the factorial design specification (one sample t
% test) and estimation and estimate Contrats.

% Mean Struct, MeanMask, Factorial design specification and
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
data.outputDir = [pathIn '\analyses\fmri\Analyses\group\task\rfx_' date '\' name extension '\'];
data.anatDir = '\anat\';
mkdir([ pathIn '\analyses\fmri\Analyses\group\task\rfx_' date '\' name extension])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matlabbatch = {};
% Create Mean Structural Image

fprintf(1,'BUILDING JOB: Create Mean Structural Image...')
for idx_sub = 1:size(data.subjects,2)
    StructFile = deblank({(ls([data.dir char(data.subjects{idx_sub}) '\task\' data.anatDir 'w*.nii']))});
    matlabbatch{1}.spm.util.imcalc.input{idx_sub,:} = [data.dir char(data.subjects{idx_sub})  '\task\' data.anatDir '\' StructFile{1}];

end

matlabbatch{1}.spm.util.imcalc.output = 'MeanStruct.img';
matlabbatch{1}.spm.util.imcalc.outdir = cellstr( data.outputDir);
matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3+i4+i5+i6+i7+i8+i9+i10+i11+i12+i13+i14+i15+i16+i17+i18+i19+i20+i21+i22+i23+i24+i25+i26+i27+i28)/28'; % NUMBER OF SUBJECTS INCLUDED
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;


% Create Mask Image
fprintf(1,'BUILDING JOB: Create Mask Image...')

for idx_sub = 1:size(data.subjects,2)
    MaskFile = deblank({(ls([data.dir char(data.subjects{idx_sub}) '\task\ffx\analysisFfx_' date '\' name extension '\mask.nii']))});
    matlabbatch{2}.spm.util.imcalc.input{idx_sub,:} = [data.dir char(data.subjects{idx_sub}) '\task\ffx\analysisFfx_' date '\' name extension '\' MaskFile{1}];
end

matlabbatch{2}.spm.util.imcalc.output = 'MeanMask.nii';
matlabbatch{2}.spm.util.imcalc.outdir = cellstr(data.outputDir);
matlabbatch{2}.spm.util.imcalc.expression = '(i1+i2+i3+i4+i5+i6+i7+i8+i9+i10+i11+i12+i13+i14+i15+i16+i17+i18+i19+i20+i21+i22+i23+i24+i25+i26+i27+i28)>0.75*28'; % NUMBER OF SUBJECTS INCLUDED

matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{2}.spm.util.imcalc.options.mask = 0;
matlabbatch{2}.spm.util.imcalc.options.interp = 1;
matlabbatch{2}.spm.util.imcalc.options.dtype = 4;


eval (['save ' data.dir '\JOBS\jobs_Means_SPM12_CLTMR_MNI' 'matlabbatch'])
fprintf(1,'Create Mean Struct and Mask IMAGES...')
spm_jobman('run',matlabbatch)

matlabbatch = {};
% Factorial design specification
fprintf(1,'BUILDING JOB: Factorial Design Specification')

%Load la liste des contrastes d'interet par analyse RFX
if isempty(extension)
    load([pathIn '\analyses\fmri\Analyses\sub-01\task\ffx\analysisFfx_' date '\contrast\contrast_' name '.mat'])
else
    load([pathIn '\analyses\fmri\Analyses\sub-01\task\ffx\analysisFfx_' date '\contrast\PPI_analysis\contrastes.mat'])
end

%RFX CONTRAST LIST
for idx_contrast = 1 : size(contrastes,2)
    listContrast= {};
    for idx_sub = 1 : size(data.subjects,2)%only one group in PoC!
        ConList = dir([data.dir char(data.subjects{idx_sub})  '\task\ffx\analysisFfx_' date '\' name '\' sprintf('scon_%0.4d.nii',idx_contrast)]);
        listContrast.path{idx_sub,:} = [data.dir char(data.subjects{idx_sub})  '\task\ffx\analysisFfx_' date '\' name '\' ConList.name];
        matlabbatch{idx_contrast}.spm.stats.factorial_design.des.t1.scans{idx_sub,:} = [data.dir char(data.subjects{idx_sub})  '\task\ffx\analysisFfx_' date '\' name extension '\' ConList.name]; % t1: One sample T Test - t2  Two sample T Test %MG changed to t1

    end
    
    matlabbatch{idx_contrast}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{idx_contrast}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{idx_contrast}.spm.stats.factorial_design.masking.em = cellstr(fullfile( data.outputDir, 'MeanMask.nii'));
    matlabbatch{idx_contrast}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{idx_contrast}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{idx_contrast}.spm.stats.factorial_design.globalm.glonorm = 1;
    
    
    matlabbatch{idx_contrast}.spm.stats.factorial_design.dir =  cellstr([data.outputDir '\' sprintf('scon_%0.4d.nii',idx_contrast)]);


end

eval (['save ' data.dir '/JOBS/jobs_FDS_SPM12_CLTMR_MNI' 'matlabbatch'])
fprintf(1,'Factorial Design Specification...')
spm_jobman('run',matlabbatch)
matlabbatch = {};

% Factorial design estimation
fprintf(1,'BUILDING JOB: Factorial Design Estimation')

for idx_contrast = 1 : size(contrastes,2)
    
    %%%%
    matlabbatch{idx_contrast}.spm.stats.fmri_est.spmmat = cellstr([data.outputDir '\' sprintf('scon_%0.4d.nii',idx_contrast) '\SPM.mat']); 
    matlabbatch{idx_contrast}.spm.stats.fmri_est.method.Classical = 1;
end

eval (['save ' data.dir '/JOBS/jobs_FDE_SPM12_CLTMR_MNI' 'matlabbatch'])
fprintf(1,'Factorial Design Estimation...')
spm_jobman('run',matlabbatch)
matlabbatch = {};

%Contrast estimation
fprintf(1,'BUILDING JOB: Contrast estimation')

for idx_contrast = 1 : size(contrastes,2)
    
    matlabbatch{idx_contrast}.spm.stats.con.spmmat = cellstr([data.outputDir '\' sprintf('scon_%0.4d.nii',idx_contrast) '\SPM.mat']); 
    matlabbatch{idx_contrast}.spm.stats.con.consess{1}.tcon.name = '+T';
    matlabbatch{idx_contrast}.spm.stats.con.consess{1}.tcon.convec = [1];
    matlabbatch{idx_contrast}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{idx_contrast}.spm.stats.con.consess{2}.tcon.name = '-T';
    matlabbatch{idx_contrast}.spm.stats.con.consess{2}.tcon.convec = [-1];
    matlabbatch{idx_contrast}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{idx_contrast}.spm.stats.con.delete = 0;
    
end

eval (['save ' data.dir '/JOBS/jobs_Contrasts_SPM12_CLTMR_MNI' 'matlabbatch'])
fprintf(1,'Contrast Estimation...')
spm_jobman('run',matlabbatch)
matlabbatch = {};



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%