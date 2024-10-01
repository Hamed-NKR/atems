clc
clear
close all
warning('off')

%% load previously saved image variables

fname = '28AUG24-HAL-End-Slider'; % name of the MATLAB worksapce file that has aggregate info
fdir = 'D:\HN\AUG24Onward\TEM\New3\ATEMS_Area'; % directory to the file to be imported

lbl = 'PFA_ET_NIT_28AUG24_HAL_End'; % label to be added later to the end of variables

fadd = cell2mat(strcat(fdir, {'\'}, fname,'.mat')); % load the MATLAB workspace file

load(fadd);

% ii = sort(randperm(length(Aggs), 25)); % aggregate ids to be looped over for labeling (empty means all)

%% assign points and id labels to the aggregates in images

n_agg = length(Aggs); % number of aggregates
randpnt = cell(n_agg,1); % initialize placeholder for random sampling points
fname_out = cell(n_agg,1); % address to the output images

if ~exist('ii', 'var') || isempty(ii)
    ii = 1 : n_agg;
end

k = 1;

while k <= length(ii) % loop over selected aggregates that are already segmented from images
    
    i = ii(k); 

    % Initialize the satisfaction variable
    is_satisfied = false;

    % Loop until the user is satisfied
    while ~is_satisfied

        % load images with aggregates being highlighted
        f1 = figure;
        j = i;
        while isempty(Aggs(j).image)
            j = j - 1;
        end
        tools.imshow_binary(Aggs(j).image, Aggs(i).binary);

        % get user's input for number of reandom points to be sampled
        n_samp = input('Please enter an integer for the number of random sampling points: ');

        % identify and randomize points within a binary image
        randpnt{i} = find(Aggs(i).binary);
        randpnt{i} = randpnt{i}(randperm(length(randpnt{i})));

        % sample n_samp random points within the binary image
        randpnt{i} = datasample(randpnt{i}, n_samp, 'Replace', false);

        % convert the point indices to row and vector locations
        [xpos1, ypos1] = ind2sub(size(Aggs(i).binary), randpnt{i});
        randpnt{i} = [xpos1, ypos1];
        
        f2 = figure;
        tools.imshow(Imgs(Aggs(i).img_id).raw);
        truesize;  % Ensure the image is displayed at true size, without interpolation
        hold on

        % plot sampling points on top of binary images
        plot(ypos1, xpos1, 'r.', 'MarkerSize', 1);

        [xpos2, ypos2] = ginput(1); % asssign points for printing aggregate id number

        text(xpos2, ypos2, num2str(Aggs(i).id), 'Color', 'red',...
        'FontSize', 14, 'FontWeight', 'bold'); % print aggregate id

        % Prompt user for feedback
        response = input('Are you satisfied with the result? (y/n): ', 's');

        if strcmpi(response, 'y')
            is_satisfied = true;
        else 
            close(f2)
        end

    end
    
    % file name for output image
    Aggs(i).fname_out = strcat(lbl, '_', num2str(Aggs(i).id));
    
    % save image array, image filename, and number of sampling points in...
    %   ...the aggregate structure
    Aggs(i).image_labeled = frame2im(getframe(gcf));
    Aggs(i).n_samp = n_samp;

    if ~exist(cell2mat(strcat('morphout', {'\'}, lbl)), 'dir')
        mkdir(cell2mat(strcat('morphout', {'\'}, lbl)));
    end
    
    fprintf('Exporting...\n\n')
    exportgraphics(f2, cell2mat(strcat('morphout', {'\'}, lbl, {'\'},...
        Aggs(i).fname_out, '.png')), 'BackgroundColor','none',...
        'ContentType','vector', 'Resolution', 2500)
    
    close all

    k = k + 1;
    
end

%% add test label  to the end of variables to be exported

% clear redundant variables
clear i j k xpos1 ypos1 xpos2 ypos2 n_samp n_agg response is_satisfied

vars = who; % call variable names
vars(strcmp(vars,'lbl')) = []; % exclude the labeling variable

for i = 1 : length(vars)

    varname = vars{i};  % Get the name of the variable as a string

    % Modify the variable name by adding a suffix
    newVarName = [varname '_' lbl];

    % Assign the value of the original variable to the new variable
    eval([newVarName ' = ' varname ';']);  % Dynamically create a new variable

    % Delete the old variable
    clear(varname);

end

clear i vars varname newVarName

save(cell2mat(strcat('morphout', {'\'}, lbl, {'\'}, lbl, '_vars.mat')))




