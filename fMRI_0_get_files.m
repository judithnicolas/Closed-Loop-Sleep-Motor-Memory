function o_files = fMRI_0_get_files(i_folder, i_ext, i_prefix)
%
%   o_files = get_files(i_folder, ext)
%   
%   i_folder: [string] Folder to list
%   i_ext:    [string] Extension file
%   i_prefix: [string] Begining of the filename
%   o_files:  [char]   List of files found
% 
% abore: 10 septembre 2015
%   - creation get_files
% abore: 10 septembre 2015
%   - modification of input file - i_
%   - add documentation
% abore: 17 septembre 2015
%   - add i_prefix
%

if nargin < 3 i_prefix=''; end;

[o_files] = spm_select('List', i_folder, ['^', i_prefix , '.*' , i_ext ,'$']);

if isempty(o_files)
    error('MyComponent:missingObject',...
    ['No files found in ', i_folder, ' folder']);
else
    o_files = strcat(i_folder, o_files);
end

