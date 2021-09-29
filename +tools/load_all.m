% LOAD_ALL          Navigate to a folder and load all mat files.
% 
% INPUTS: NONE.
%
% OUTPUTS: 
% [FOLDERNAME] = Path containing the folder name and the 2 folders above.
% [H] = Output struct containing names and data of each mat file.
%
%  Darwin Zhu, 2021-06-09
%=========================================================================%

function [foldername, h] = load_all()

% Preserve the original location
start_path = pwd;

% Start UI, navigate to folder
path = uigetdir;

% Get folder name, and the 2 directories above it.
[folderpartsrest, folderparts1] = fileparts(path);
[folderpartsrest, folderparts2] = fileparts(folderpartsrest);
[folderpartsrest, folderparts3] = fileparts(folderpartsrest);
foldername = [folderparts3,'/', folderparts2,'/', folderparts1];
cd(path);
files = dir('*.mat');

a = size(files);
h = struct([]);

% Iterate over files struct
for i = 1:a(1)
    val = load(files(i).name);
    fields = fieldnames(val);
    
    %Store the name and struct file into the h struct.
    h(i).name = fields{1,1};
    h(i).value = val.(fields{1,1});
end

% Navigate back to original folder.
cd(start_path);

end