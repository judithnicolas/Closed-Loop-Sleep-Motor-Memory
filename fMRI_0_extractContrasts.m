function  [C,contrastes] = fMRI_14_extractContrasts(pathIn,sub,date,mod)%% Import data from text file

%%

load([pathIn '\analyses\fmri\Analyses\' sub '\task\ffx\analysisFfx_' date '.mat\SPM.mat'])
load([pathIn '\analyses\fmri\scripts\contrastList.mat'])

%%
if mod ==0
    cond1 = cellfun(@isempty,regexp(SPM.xX.name,'Median')); % does NOT contain Median
    name = 'raw';
else
    cond1 = ~cellfun(@isempty,regexp(SPM.xX.name,'Median')); % does contain Median
    name = 'mod';
end

cond2 = cellfun(@isempty,regexp(SPM.xX.name,' R')); % do not contain R
cond3 = cellfun(@isempty,regexp(SPM.xX.name,' bin')); % do not contain bin
cond4 = cellfun(@isempty,regexp(SPM.xX.name,' constant')); % do not contain bin

condAll = cond1.*cond2.*cond3.*cond4;

%%
C = zeros(size(contrastList,1)*2,size(SPM.xX.X,2)); % TOTAL NUMBER OF CONTRASTS TO RUN FOR MAIN EFFECTS 
contrastes = struct('C',[],'name',[]);

for idx_contrasts = 1 : size(contrastList,1)
    
    contrastName = contrastList(idx_contrasts,:);
    [C,contrastes] =  extractValuesContrast(SPM,C,contrastes,condAll,contrastName);

end
mkdir([pathIn '\analyses\fmri\Analyses\' sub '\task\ffx\analysisFfx_' date '.mat\contrast\'])
save([pathIn '\analyses\fmri\Analyses\' sub '\task\ffx\analysisFfx_' date '.mat\contrast\contrast_' name '.mat'],'C' , 'contrastes')

