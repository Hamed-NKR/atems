function [Aggs, imgs_binary] = categorize_manu(imgs, pixsizes, opts, fname)

% CATEGORIZE_MANU manually segments the aggregates and then asks the...
%   ...user for inputs on morphological types.
% ----------------------------------------------------------------------- %
% imgs: a cell array of "cropped" grayscale images to be analyzed
% pixsizes: image pixel sizes
% opts: a structure of the function's options
% fname: file-names of the images
% Aggs: particles' table of properties
% imgs_binary: particles' segmented binary images
% ----------------------------------------------------------------------- %

%% Part1: Categorize the particles %%

n_imgs = length(imgs); % number of images

% initialize inputs if not given
if ~(exist('opts', 'var') && isfield(opts, 'sizing') &&...
        isfield(opts, 'ui'))
    opts = struct('sizing', [], 'ui', []);
end

opts_sizing = opts.sizing;
opts_ui = opts.ui;

if isempty(opts_ui)
    opts_ui = 'on'; % default to get the user's morphological inputs
end

if isempty(opts_sizing)
    opts_sizing = 'manu'; % default to be manual aggregate sizing
end

if ~exist('fname','var'); fname = []; end
if isempty(fname); fname = cell(length(imgs)); end

% manually segment the aggregates
if ismember(opts_sizing, {'MANU', 'Manu', 'manu'})
    imgs_binary = agg.seg_slider(imgs);
elseif ismember(opts_sizing, {'AUTO', 'Auto', 'auto'})
    imgs_binary = agg.seg_kmeans(imgs, pixsizes);
else
    error(['Invalid sizing option input!', newline,...
        'Should be ''Manu''/''Auto''.'])
end

% obtain aggregate properties
Aggs = agg.analyze_binary(imgs_binary, pixsizes, imgs, fname);

agg_id = cat(1, Aggs.id); % aggregate ids
img_id = cat(1, Aggs.img_id); % image ids corresponding to the aggregates
n_agg_i = zeros(n_imgs,1); % initialize number of aggregates within each...
    % ...image

if ismember(opts_ui, {'ON', 'On', 'on'})
    opts_ui_bin = 1; % gui checking handle
    
    % initialize the figure properties for TEM images to be plotted 
    figure;
    hf0 = gcf;
    hf0.Position = [0, 0, 800 * 2, 800]; % Position and size
    set(hf0, 'color', 'white'); % Background color
    
    % set the figure layout 
    tt0 = tiledlayout(1, 2);
    tt0.TileSpacing = 'compact';
    tt0.Padding = 'compact';
    
elseif ismember(opts_ui, {'OFF', 'Off', 'off'})
    opts_ui_bin = 0;
    morphstr = {'fs', 'cs', 'tb', 'sb', 'h', 'm'};
    
else
    error(['Invalid user interface option input!', newline,...
        'Should be ''On''/''Off''.'])
end

for i = 1 : n_imgs
    
    ii = sum(n_agg_i(1 : i-1)) + 1; % global id of the first aggregate...
        % ...within the image
    n_agg_i(i) = length(agg_id(img_id == i));
    
    if opts_ui_bin
        tt0_a = nexttile;
        tools.imshow(Aggs(ii).image); % display the original image
        title(tt0_a, 'Original TEM image')
        
        tt0_b = nexttile;
    end
    
    for j = 1 : n_agg_i(i)
        jj = sum(n_agg_i(1 : i-1)) + j; % aggregate's global id

        if opts_ui_bin
            % shade the isolated aggregate area
            cmap = ones(max(max(Aggs(jj).binary)), 1) * [0.12, 0.59, 0.96];
            im2 = labeloverlay(Aggs(ii).image, Aggs(jj).binary,...
                'Transparency', 0.6, 'Colormap', cmap);
        
            % display the aggregate edge
            agg_edge = edge(Aggs(jj).binary, 'sobel'); % capture the edge
            se = strel('disk',1);
            edge_dilated = imdilate(agg_edge, se); % strengthen the outline
            im2 = uint8(~edge_dilated) .* im2; % add borders to labeled regions

            tools.imshow(im2); % display the highlighted aggregate and the edge
            title(tt0_b, 'Highlighted image')
        end
        
        morph_check = 0; % initialize the checking variable for...
            % ...the user's inputs
        
        while ~all(morph_check) % check if acceptable inputs are provided
            
            if opts_ui_bin
                % request inputs on circulariry and optical depth from the user
                prompt = {['\fontsize{12}Please enter the',...
                    '\fontsize{12}\bfMorphological Type:', newline,...
                    '\fontsize{12}\rm\it(FS/Fs/fs: Fractal soot, Compact soot,',...
                    '\fontsize{12}\rm\it TB/Tb/tb: Tarball, SB/Sb/sb: Softball, H/h: Hybrid, M/m: Miscellaneous)']};
                dlgtitle = 'User''s morphological inputs';
                dims = [1 100];
                defaultans = {'', '', '', ''};
                dlgopts.Interpreter = 'tex';
                morph = inputdlg(prompt, dlgtitle, dims, defaultans, dlgopts);
                
            else
                morph = morphstr{randsample(6,1)};
            end
            
            % feed the inputs to the aggregate properties cell array
            switch morph % morph. type
                case {'FS', 'Fs', 'fs'}
                    Aggs(jj).Type = 'Fractal soot';
                case {'CS', 'Cs', 'cs'}
                    Aggs(jj).Type = 'Compact soot';
                case {'TB', 'Tb', 'tb'}
                    Aggs(jj).Type = 'Tarball';
                case {'SB', 'Sb', 'sb'}
                    Aggs(jj).Type = 'Softball';
                case {'H', 'h'}
                    Aggs(jj).Type = 'Hybrid';
                case {'M', 'm'}
                    Aggs(jj).Type = 'Miscellaneous';
                otherwise
                    % issue warning
                    warn = questdlg(['Invalid input for morph. type!',...
                        newline, '(please see the descriptions in the prompt input window)'],...
                        'warning', 'Continue', 'Cancel', 'Continue');
                    % handle response
                    switch warn
                        case 'continue'
                            continue
                        case 'cancel'
                            warning(['Invalid input for morph. type!', newline,...
                                '(please see the descriptions in the prompt input window)'])
                            return
                    end
            end
            morph_check = 1; % morph. type input is valid
            
            if opts_ui_bin
                % clear the highlighted plot to prepare for next aggregate
                ha0 = gca;
                cla(ha0)
            end
        end
    end
    
    if opts_ui_bin
        clf(hf0) % clear figure to prepare for the next image
    end
end

if opts_ui_bin
    close(hf0) % close figure
    clear hf0 % clear figure handle
end
end
