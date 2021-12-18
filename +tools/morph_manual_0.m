function [ff, hf] = morph_manual_0(Aggs)
% MORPH_MANUAL characterizes the aggregates based on their morphological...
%   ...properties including circularity and optical depth
% ----------------------------------------------------------------------- %
% Aggs: aggregate data structure
% ff: morphological frequencies
% hf: results figure handle
% ----------------------------------------------------------------------- %

agg_id = cat(1, Aggs.id); % aggregate ids
img_id = cat(1, Aggs.img_id); % image ids
n_img = length(unique(img_id)); % number of images
n_agg_i = zeros(n_img,1); % initialize number of aggregates within each...
    % ...image
n_agg_tot = length(Aggs); % total number of aggregates across all images

% initialize the figure properties for TEM images to be plotted 
figure;
hf0 = gcf;
hf0.Position = [0, 0, 900 * 2, 900]; % Position and size
set(hf0, 'color', 'white'); % Background color

% set the figure layout 
tt0 = tiledlayout(1, 2);
tt0.TileSpacing = 'compact';
tt0.Padding = 'compact';

ff0 = zeros(n_agg_tot, 3, 2); % initialize the fuzzy frequency array;... 
    % ...rows: image repository; columns: morphological level (H/M/L);...
    % ...layers: circularity and optical depth)

for i = 1 : n_img
    ii = sum(n_agg_i(1 : i-1)) + 1; % global id of the first aggregate...
        % ...within the image
    n_agg_i(i) = length(agg_id(img_id == i));
    
    tt0_a = nexttile;
    tools.imshow(Aggs(ii).image); % display the original image
    title(tt0_a, 'Original TEM image')
    
    tt0_b = nexttile;
    for j = 1 : n_agg_i(i)
        jj = sum(n_agg_i(1 : i-1)) + j; % aggregate's global id

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
        
        morph_check = [0, 0]; % initialize the checking variable for...
            % ...the user's inputs
        
        while ~all(morph_check) % check if acceptable inputs are provided
            
            % request inputs on circulariry and optical depth from the user
            prompt = {['Enter circularity level: ', newline,...
                '(H/h: High [circle], M/m: Moderate, L/l: Low [irregular]):'],...
                ['Enter optical depth level: ', newline,...
                '(H/h: High [dark], M/m: Moderate, L/l: Low [bright]):']};
            dlgtitle = 'User''s morphological inputs';
            dims = [1 60];
            morph = inputdlg(prompt, dlgtitle, dims);
            
            % feed the inputs to the aggregate properties structure
            switch morph{1}
                case {'L', 'l'}
                    ff0(jj,1,1) = 1; % relatively acircular
                case {'M', 'm'}
                    ff0(jj,2,1) = 1; % moderately acircular
                case {'H', 'h'}
                    ff0(jj,3,1) = 1; % aggregate to be nearly circular
                otherwise
                    % issue warning
                    warn1 = questdlg(['Invalid input for circularity!',...
                        newline, '(please see the descriptions in the prompt input window)'],...
                        'warning', 'continue', 'cancel', 'continue');
                    % Handle response
                    switch warn1
                        case 'continue'
                            continue
                        case 'cancel'
                            warning(['Invalid input for circularity!', newline,...
                                '(please see the descriptions in the prompt input window)'])
                            return
                    end
            end
            morph_check(1) = 1; % circularity input is valid
            
            switch morph{2}
                case {'L', 'l'}
                    ff0(jj,1,2) = 1; % relatively bright
                case {'M', 'm'}
                    ff0(jj,2,2) = 1; % moderately dark
                case {'H', 'h'}
                    ff0(jj,3,2) = 1; % aggregate to be nearly black
                otherwise
                    % issue warning
                    warn2 = questdlg(['Invalid input for circularity!',...
                        newline, '(please see the descriptions in the prompt input window)'],...
                        'warning', 'continue', 'cancel', 'continue');
                    % Handle response
                    switch warn2
                        case 'continue'
                            continue
                        case 'cancel'
                            warning(['Invalid input for optical depth!', newline,...
                                '(please see the descriptions in the prompt input window)'])
                            return
                    end
            end
            morph_check(2) = 1; % optical depth input is valid
            
            % clear the highlighte plot to prepare for next aggregate
            ha0 = gca;
            cla(ha0)
        end
    end
    
    clf(hf0) % clear figure to prepare for the next image
end

ff = [nnz(ff0(:,1,1)), nnz(ff0(:,2,1)),nnz(ff0(:,3,1));...
    nnz(ff0(:,1,2)), nnz(ff0(:,2,2)), nnz(ff0(:,3,2))]; 
        % total (image integrated) frequencies (row1: circulariry,...
            % ...row2: optical depth)
ff = ff * 100 / n_agg_tot; % convert to percentage

close(hf0) % close figure
clear hf0 % clear figure handle

% initialize the figure for the manual morphological results
figure;
hf = gcf;
hf.Position = [0, 0, 700 * 2, 700]; % Position and size
set(hf, 'color', 'white'); % Background color

% set the figure layout 
tt = tiledlayout(1, 2);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

% Define colormaps
cm1 = colormap(autumn);
cm2 = colormap(gray);

% circularity frequencies
nexttile;
x1 = categorical({'Irregular', 'Moderate', 'Circular'});
x1 = reordercats(x1,{'Irregular', 'Moderate', 'Circular'});
bc1 = bar(x1, ff(1,:), 'FaceColor', 'flat');
bc1.CData = [cm1(200,:); cm1(110,:); cm1(20,:)];
box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 12)
ylim([0 100])
xlabel('c (-)', 'FontName', 'SansSerif', 'FontSize', 14,...
    'FontWeight', 'bold')
ylabel('f_n (%)', 'FontName', 'SansSerif', 'FontSize', 14,...
    'FontWeight', 'bold')

% optical depth frequencies
nexttile;
x2 = categorical({'Light', 'Moderate', 'Dark'});
x2 = reordercats(x2,{'Light', 'Moderate', 'Dark'});
bc2 = bar(x2, ff(2,:), 'FaceColor', 'flat');
bc2.CData = [cm2(200,:); cm2(110,:); cm2(20,:)];
box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 12)
ylim([0 100])
xlabel('z_o_p_t (-)', 'FontName', 'SansSerif', 'FontSize', 14,...
    'FontWeight', 'bold')
ylabel('f_n (%)', 'FontName', 'SansSerif', 'FontSize', 14,...
    'FontWeight', 'bold')

title(tt, 'Visually interpreted morphological properties',...
    'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 16)

if nargout < 2
    clear hf % delete the results figure handle if not requested as an...
        % ...output
end

end

