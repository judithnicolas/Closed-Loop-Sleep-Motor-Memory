function o_matlabbatch = fMRI_0_coregister(anat, mean, other)
%
%   function o_matlabbatch = coregister(i_ref, i_source)
%
%   i_ref:          [string]    Reference for coregistration  (steady)
%   i_source:       [string]    Source for coregistration
%
%   abore: 17 Septembre 2015
%       - creation of coregister
o_matlabbatch = [];

o_matlabbatch{end+1}.spm.spatial.coreg.estimate.ref = cellstr(anat);
o_matlabbatch{end}.spm.spatial.coreg.estimate.source(1) = cellstr(mean);
o_matlabbatch{end}.spm.spatial.coreg.estimate.other = cellstr(other);
%     o_matlabbatch{end}.spm.spatial.coreg.estwrite.other = cfg_dep('Expand image frames: Expanded filename list.', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
o_matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
o_matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
o_matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 ...
    0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
o_matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
