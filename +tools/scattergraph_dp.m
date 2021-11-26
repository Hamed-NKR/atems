
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

% offsets on text points
dx = 0.1;
dy = 0.1;

% initialize fields
dp = [];
dp_manual = [];
dp_std = [];
err = struct([]);
k = 1;

ids = [];

% Iterate over struct h and compile values into dp, dp_manual, and dp_std.
for i = 1:max(len(1),len(2))
    h_val = h(i).value;
    h_val_size = size(h_val);
    h_val_name = h(i).name;
    % Iterate over the h_val element in h. 
    for j = 1:h_val_size(2)
        dp(end+1) = h_val(j).dp;
        dp_manual(end+1) = h_val(j).dp_manual;
        dp_std(end+1) = h_val(j).dp_std;
        ids(end+1) = h_val(j).id;
        % If value of dp exceeds 50% of expected, store its name & id in
        % err array.
        if (dp_manual(end) > dp(end)*1.5 || dp_manual(end) < dp(end)*0.5)
            err(k).name = h_val_name;
            err(k).id = h_val(j).id;
            err(k).dp = h_val(j).dp;
            err(k).dp_manual = h_val(j).dp_manual;
            k = k + 1;
        end
    end
end

% Create an errorbar plot. 
errorbar(dp, dp_manual, dp_std,'vertical', 'ko', 'LineWidth', 2);
xlabel('PCM dp (nm)');
ylabel('Manual median dp (nm)');
%title({append(foldername, ' - Manual vs PCM dp')}, 'Interpreter', 'none');
title ({'Comparing manual to PCM dp'});

hold on;

% Plot 1:1 expected line in red
maxaxis = max(max(dp), max(dp_manual));
plot([0,maxaxis],[0,maxaxis], 'r-', 'LineWidth', 2);

% Plot best fit line in blue
p = polyfit(dp, dp_manual,1);
y1 = polyval(p, dp);
plot(dp, y1, 'b-', 'LineWidth', 2);
set(gca, 'FontSize', 20);

legend({'Particle dp', '1:1 ratio', 'Polyfit best fit line'}, 'Location', 'southeast')

% Plot id numbers onto points

% for k = 1:length(ids)
%     text(dp(k) + dx, dp_manual(k) + dy, string(ids(k)));
% end
hold off;

end