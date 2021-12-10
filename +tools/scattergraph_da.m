
% SCATTERGRAPH      Draws a scatter plot comparing dp and dp_manual of all
%                   .mat files in a folder, with error bars. Outputs a
%                   struct containing the names and ids of aggregates with
%                   >50% error. Prints the ids of >50% error aggregates
%                   onto scatter plot. 
%
% INPUTS: NONE
%
% OUTPUTS: [ERR] = Struct containing names and ids of all primary particles
% that fall outside of the 50% expected range between dp_manual and pcm.
% Author:           Darwin Zhu, 2021-06-09

%=========================================================================%

function [err] = scattergraph_dp()
% Run load_all to compile all case files into struct h.
[foldername, h] = tools.load_all();
len = size(h);
% Load kmeans file using ui
    uiwait(msgbox('Please navigate to the .mat file containing the kmeans struct',...
            'Loading aggregate struct...','help'));
    [kmeans_file, kmeans_folder] = uigetfile('*.mat');
    kmeans_full = fullfile(kmeans_folder, kmeans_file);
    
    
    % Check if file was correctly selected.
    if exist(kmeans_full, 'file')
        % Normal situation - they picked an existing file.
        kmeans_in = load(kmeans_full);
        kmeans_aggname = fieldnames(kmeans_in);
        Aggs_kmeans = kmeans_in.(kmeans_aggname{1});
        % Now do something with storedStructure, like extract fields into new variables or whatever you want.
    else
        % Error: Would only get here if they typed in a name of a non-existant file
        % instead of picking one from the folder.
        warningMessage = sprintf('Warning: mat file does not exist:\n%s', kmeans_full);
        uiwait(errordlg(warningMessage));
        return;
    end
len_kmeans = size(Aggs_kmeans);
% offsets on text points
dx = 0.1;
dy = 0.1;

% initialize fields
da = [];
da_kmeans = [];
err = struct([]);
k = 1;

ids = [];

% Iterate over struct h and compile values into dp, dp_manual, and dp_std.
for i = 1:max(len(1),len(2))
    h_val = h(i).value;
    h_val_size = size(h_val);
    h_val_name = h(i).name;
    if (i <= len_kmeans(2))
        da_kmeans(end+1) = Aggs_kmeans(i).da;
    end
    % Iterate over the h_val element in h. 
    for j = 1:h_val_size(2)
        da(end+1) = h_val(j).da;
        ids(end+1) = h_val(j).id;

    end
end

% Create an errorbar plot. 
scatter(da(1:len_kmeans(2)), da_kmeans, 'rx', 'LineWidth', 2);
hold on;

ylabel('kmeans da (nm)');
xlabel('Manual da (nm)');
%title({append(foldername, ' - Manual vs PCM dp')}, 'Interpreter', 'none');
title ({'Comparing manual to PCM da'});


% Plot 1:1 expected line in red
maxaxis = max(max(da), max(da_kmeans));
plot([0,maxaxis],[0,maxaxis], 'r-', 'LineWidth', 2);

% Plot best fit line in blue
%p = polyfit(dp, dp_manual,1);
%y1 = polyval(p, dp);
%plot(dp, y1, 'b-', 'LineWidth', 2);
set(gca, 'FontSize', 20);

legend({'kmeans seg da', '1:1 ratio'}, 'Location', 'southeast')

% Plot id numbers onto points

% for k = 1:length(ids)
%     text(dp(k) + dx, dp_manual(k) + dy, string(ids(k)));
% end
hold off;

end