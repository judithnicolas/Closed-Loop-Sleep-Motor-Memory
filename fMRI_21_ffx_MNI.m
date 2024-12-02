function fMRI_21_ffx_MNI(sub,pathIn,date,mod)

if mod == 0 
    name = 'raw';
else
    name = 'mod';
end
load ([pathIn '\analyses\fmri\Analyses\' sub '\task\ffx\analysisFfx_' date '\contrast\contrast_' name '.mat'])

strdir = [pathIn '\analyses\fmri\Analyses\' sub '\task\ffx\analysisFfx_' date '\'];

fprintf(1,'PROCESSING SUBJECT : %s\n',sub)
mkdir(fullfile(strdir,'\',name))
copyfile (fullfile(strdir,'SPM.mat'),fullfile(strdir,'\',name,'\'))
copyfile (fullfile(strdir,'\beta*'),fullfile(strdir,'\',name,'\'))
copyfile (fullfile(strdir,'\ResMS*'),fullfile(strdir,'\',name,'\'))
copyfile (fullfile(strdir,'\mask*'),fullfile(strdir,'\',name,'\'))
copyfile (fullfile(strdir,'\beta*'),fullfile(strdir,'\',name,'\'))
matlabbatch{3} = {};
matlabbatch{3}.spm.stats.con.spmmat = cellstr(fullfile(strdir,'\',name,'\','SPM.mat'));

for icon = 1:size(contrastes,2)
    matlabbatch{3}.spm.stats.con.consess{icon}.tcon.name = contrastes(icon).name;
    matlabbatch{3}.spm.stats.con.consess{icon}.tcon.convec = contrastes(icon).C;
    matlabbatch{3}.spm.stats.con.consess{icon}.tcon.sessrep = 'none';
end % icon

matlabbatch{3}.spm.stats.con.delete = 0;

eval (['save ' pathIn '\analyses\fmri\Analyses\JOBS\jobs_FFX_SPM12_CLTMR_MNI_' name '_' sub '.mat ' ' matlabbatch'])

% Pruning the empty options
for ii = 1:2
    matlabbatch(1) = []; % no curly brackets to remove array
end

spm_jobman('run',matlabbatch)


clear matlabbatch
