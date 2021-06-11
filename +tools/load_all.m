
% LOAD_ALL          Navigate to a folder and load all mat files to a struct
% Author:           Darwin Zhu, 2021-06-09

%=========================================================================%

function [h] = load_all()

% Preserve the original location
start_path = pwd;

% Start UI, navigate to folder
path = uigetdir;
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