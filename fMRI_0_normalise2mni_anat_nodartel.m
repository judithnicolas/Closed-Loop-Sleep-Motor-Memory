function o_matlabbatch = fMRI_0_normalise2mni_anat_nodartel(i_flowfield, i_files2normalise)
% 
%   function o_matlabbatch = normalise2mni(i_flowfield, i_files2normalise)
% 
%   i_flowfield:    [string]
%   i_files2normalise:  [string]
% 
%   o_matlabbatch: [array]
% 
%   ga: 23 July 2018
%       - creation of normalise2mni_anat_nodartel
% 

o_matlabbatch = [];

o_matlabbatch{end+1}.spm.spatial.normalise.write.subj.def(1) = cellstr(i_flowfield);
o_matlabbatch{end}.spm.spatial.normalise.write.subj.resample = cellstr(i_files2normalise);
o_matlabbatch{end}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
o_matlabbatch{end}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
o_matlabbatch{end}.spm.spatial.normalise.write.woptions.interp = 4;
o_matlabbatch{end}.spm.spatial.normalise.write.woptions.prefix = 'w';
