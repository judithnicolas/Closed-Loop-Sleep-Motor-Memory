function fMRI_0_reorient(i_file)
% 
%   function reorient(i_file)
% 
%   i_file: [string]
% 
%   bore: 17 Septembre 2015
%       - creation of reorient
% 
% Function to automatically reorient a T2 structural image on the template
% in order further proceed with the segmentation/normalisation.
% very useful as segmentation/normalisation is rather sensitive on the
% starting orienation of the image

if nargin<1
    i_file=spm_select(inf,'image');
end
tmpl=fullfile(spm('dir'),'canonical\avg152T1.nii');
vg=spm_vol(tmpl);
flags.regtype='rigid';
for i=1:size(i_file,1)
    f=strtrim(i_file(i,:));
    spm_smooth(f,'temp.nii',[12 12 12]);
    vf=spm_vol('temp.nii');
    [M,scal] = spm_affreg(vg,vf,flags);
    M3=M(1:3,1:3);
    [u s v]=svd(M3);
    M3=u*v';
    M(1:3,1:3)=M3;
    N=nifti(f);
    N.mat=M*N.mat;
    create(N);
end

