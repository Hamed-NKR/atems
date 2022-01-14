function [ff, Aggs, imgs_binary, hf] = categorize_manu_0(imgs, pixsizes,...
    imnames, opts, fout)
% ANALYZE_MORPH characterizes the aggregates based on their morphological...
%   ...properties including circularity and optical depth
% ----------------------------------------------------------------------- %
% imgs: cell array of "cropped" grayscale images to be analyzed
% pixsizes: image pixel sizes
% imnames: image file name
% opts: options structure
% fout: output file saving address
% ff: morphological frequencies
% Aggs: aggregate properties tables
% imgs_binary: aggregate segmented binary images
% hf: results figure handle
% ----------------------------------------------------------------------- %

n_imgs = length(imgs); % number of images

if ~(exist('opts', 'var') && isfield(opts, 'sizing') &&...
        isfield(opts, 'ui'))
    opts = struct('sizing', [], 'ui', []);
end

opts_sizing = opts.sizing;
opts_ui = opts.ui;

if ~exist('opts_ui', 'var') || isempty(opts_ui)
    opts_ui = 'on'; % default to get the user's morphological inputs
end

if ~exist('opts_sizing', 'var') || isempty(opts_sizing)
    opts_sizing = 'manual'; % default to be manual aggregate sizing
end

if ~exist('fout', 'dir') || isempty(fout)
    fout = pwd; % default saving address to be in the working directory
end

% manually segment the aggregates
if ismember(opts_sizing, {'MANUAL', 'Manual', 'manual'})
    imgs_binary = agg.seg_slider(imgs);
elseif ismember(opts_sizing, {'AUTO', 'Auto', 'auto'})
    imgs_binary = agg.seg_kmeans(imgs, pixsizes);
else
    error(['Invalid sizing option input!', newline,...
        'Should be ''Manual''/''Auto''.'])
end

% obtain aggregate properties
Aggs = agg.analyze_binary(imgs_binary, pixsizes, imgs, imnames);

agg_id = cat(1, Aggs.id); % aggregate ids
img_id = cat(1, Aggs.img_id); % image ids corresponding to the aggregates
n_agg_i = zeros(n_imgs,1); % initialize number of aggregates within each...
    % ...image
n_agg_tot = length(Aggs); % total number of aggregates across all images

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
    m1 = {'h', 'm', 'l'};
    m2 = {'fs', 'ds', 'cs', 'tb', 'gb', 'hy', 'mi'};
    
else
    error(['Invalid user interface option input!', newline,...
        'Should be ''On''/''Off''.'])
end

ff0 = {zeros(n_agg_tot, 3), zeros(n_agg_tot, 3), zeros(n_agg_tot, 3),...
    zeros(n_agg_tot, 7)}; % initialize the fuzzy frequency array with...
        % ...circularity, optical depth, sharpness and type storage arrays

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
        
        morph_check = zeros(4,1); % initialize the checking variable for...
            % ...the user's inputs
        
        while ~all(morph_check) % check if acceptable inputs are provided
            
            if opts_ui_bin
                % request inputs on circulariry and optical depth from the user
                prompt = {['\fontsize{12}Please enter:', newline,...
                    '\fontsize{12}\bfCircularity level',...
                    '\fontsize{12}\rm\it (H/h: High [circle], M/m: Moderate, L/l: Low [irregular]):'],...
                    ['\fontsize{12}\bfOptical depth level',...
                    '\fontsize{12}\rm\it (H/h: High [dark], M/m: Moderate, L/l: Low [bright])'],...
                    ['\fontsize{12}\rm\bfSharpness level',...
                    '\fontsize{12}\rm\it (H/h: High [distinctive], M/m: Moderate, L/l: Low [diffused])'],...
                    ['\fontsize{12}\rm\bfMorpholgical type',...
                    '\fontsize{12}\rm\it (FS/Fs/fs: Fractal Soot, DS/Ds/ds: Dense soot, CS/Cs/cs: Coated soot,',...
                    '\fontsize{12}\rm\it TB/Tb/tb: Tarball, GB/Gb/gb: Greyball, HY/Hy/hy: Hybrid, MI/Mi/mi: Miscellaneous)']};
                dlgtitle = 'User''s morphological inputs';
                dims = [3 120];
                defaultans = {'', '', '', ''};
                dlgopts.Interpreter = 'tex';
                morph = inputdlg(prompt, dlgtitle, dims, defaultans, dlgopts);
                
            else
                morph = {m1{randsample(3,1)}, m1{randsample(3,1)},...
                    m1{randsample(3,1)}, m2{randsample(7,1)}};
            end
            
            % feed the inputs to the aggregate properties cell array
            switch morph{1} % circularity
                case {'L', 'l'}
                    ff0{1}(jj,1) = 1;
                case {'M', 'm'}
                    ff0{1}(jj,2) = 1;
                case {'H', 'h'}
                    ff0{1}(jj,3) = 1;
                otherwise
                    % issue warning
                    warn1 = questdlg(['Invalid input for circularity!',...
                        newline, '(please see the descriptions in the prompt input window)'],...
                        'warning', 'continue', 'cancel', 'continue');
                    % handle response
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
            
            switch morph{2} % optical depth
                case {'L', 'l'}
                    ff0{2}(jj,1) = 1;
                case {'M', 'm'}
                    ff0{2}(jj,2) = 1;
                case {'H', 'h'}
                    ff0{2}(jj,3) = 1;
                otherwise
                    % issue warning
                    warn2 = questdlg(['Invalid input for optical depth!',...
                        newline, '(please see the descriptions in the prompt input window)'],...
                        'warning', 'continue', 'cancel', 'continue');
                    % handle response
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
            
            switch morph{3} % sharpness
                case {'L', 'l'}
                    ff0{3}(jj,1) = 1;
                case {'M', 'm'}
                    ff0{3}(jj,2) = 1;
                case {'H', 'h'}
                    ff0{3}(jj,3) = 1;
                otherwise
                    % issue warning
                    warn3 = questdlg(['Invalid input for sharpness!',...
                        newline, '(please see the descriptions in the prompt input window)'],...
                        'warning', 'continue', 'cancel', 'continue');
                    % handle response
                    switch warn3
                        case 'continue'
                            continue
                        case 'cancel'
                            warning(['Invalid input for sharpness!', newline,...
                                '(please see the descriptions in the prompt input window)'])
                            return
                    end
            end
            morph_check(3) = 1; % sharpness input is valid
            
            switch morph{4} % morph. type
                case {'FS', 'Fs', 'fs'}
                    ff0{4}(jj,1) = 1;
                case {'DS', 'Ds', 'ds'}
                    ff0{4}(jj,2) = 1;
                case {'CS', 'Cs', 'cs'}
                    ff0{4}(jj,3) = 1;
                case {'TB', 'Tb', 'tb'}
                    ff0{4}(jj,4) = 1;
                case {'GB', 'Gb', 'gb'}
                    ff0{4}(jj,5) = 1;
                case {'HY', 'Hy', 'hy'}
                    ff0{4}(jj,6) = 1;
                case {'MI', 'Mi', 'mi'}
                    ff0{4}(jj,7) = 1;
                otherwise
                    % issue warning
                    warn4 = questdlg(['Invalid input for morph. type!',...
                        newline, '(please see the descriptions in the prompt input window)'],...
                        'warning', 'continue', 'cancel', 'continue');
                    % handle response
                    switch warn4
                        case 'continue'
                            continue
                        case 'cancel'
                            warning(['Invalid input for morph. type!', newline,...
                                '(please see the descriptions in the prompt input window)'])
                            return
                    end
            end
            morph_check(4) = 1; % morph. type input is valid
            
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

ff = {[nnz(ff0{1}(:,1)), nnz(ff0{1}(:,2)), nnz(ff0{1}(:,3,1))],...
    [nnz(ff0{2}(:,1)), nnz(ff0{2}(:,2)), nnz(ff0{2}(:,3))],...
    [nnz(ff0{3}(:,1)), nnz(ff0{3}(:,2)), nnz(ff0{3}(:,3))],...
    [nnz(ff0{4}(:,1)), nnz(ff0{4}(:,2)), nnz(ff0{4}(:,3))],...
    [nnz(ff0{4}(:,4)), nnz(ff0{4}(:,5)), nnz(ff0{4}(:,6)),...
    nnz(ff0{4}(:,7))]}; % total (image integrated) frequencies

fout = strcat(fout, '\morphout'); % set the output file folder
if ~exist('fout', 'dir')
    mkdir(fout)
end
fout = strcat(fout, '\morphdata.mat'); % set the output file name
save(fout, 'ff', 'Aggs', 'imgs_binary'); % save the results in case an...
    % ...error occurs while plotting

% initialize the figure for the manual morphological results
figure;
hf = gcf;
hf.Position = [0, 0, 890, 800]; % Position and size
set(hf, 'color', 'white'); % Background color

% set the figure layout
tt = tiledlayout(3,3);
tt.TileSpacing = 'loose';
tt.Padding = 'loose';

% define colormaps
cm1 = colormap(autumn);
cm2 = colormap(gray);
cm3 = colormap(winter);
cm4 = colormap(summer);
cm5 = colormap(spring);

ff2 = cell2mat(ff); % concatinate freqs.
err = 100 * sqrt(ff2 / n_agg_tot^2 + ff2.^2 / n_agg_tot^3) ./ ff2;
    % error propagation in percent
ff2 = 100 * ff2 / n_agg_tot; % relative freq. in percent

% circularity frequency plot
nexttile(1);
x1 = categorical({'Irregular', 'Moderate', 'Circular'});
x1 = reordercats(x1,{'Irregular', 'Moderate', 'Circular'});
bc1 = bar(x1, ff2(1:3), 'FaceColor', 'flat');
bc1.CData = [cm1(200,:); cm1(110,:); cm1(20,:)];
hold on
eb1 = errorbar(x1, ff2(1:3), err(1:3), '.');
eb1.Color = [0.1 0.1 0.1];
eb1.CapSize = 5;
eb1.MarkerSize = 0.1;
box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 12,...
    'TickLength', [0.03 0.03])

ylim([0 min(1.1 * (max(ff2(1:3) + err(1:3))), 100)])
xlabel('Circularity level', 'FontName', 'SansSerif', 'FontSize', 14,...
    'FontWeight', 'bold')
ylabel('Frequency (%)', 'FontName', 'SansSerif', 'FontSize', 14,...
    'FontWeight', 'bold')
hold off

% optical depth frequency plot
nexttile(4);
x2 = categorical({'Light', 'Moderate', 'Dark'});
x2 = reordercats(x2,{'Light', 'Moderate', 'Dark'});
bc2 = bar(x2, ff2(4:6), 'FaceColor', 'flat');
bc2.CData = [cm2(210,:); cm2(130,:); cm2(50,:)];
hold on
eb2 = errorbar(x2, ff2(4:6), err(4:6), '.');
eb2.Color = [0.6350 0.0780 0.1840];
eb2.CapSize = 5;
eb2.MarkerSize = 0.1;
box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 12,...
    'TickLength', [0.03 0.03])
ylim([0 min(1.1 * (max(ff2(4:6) + err(4:6))), 100)])
xlabel('Optical depth level', 'FontName', 'SansSerif', 'FontSize', 14,...
    'FontWeight', 'bold')
ylabel('Frequency (%)', 'FontName', 'SansSerif', 'FontSize', 14,...
    'FontWeight', 'bold')
hold off

% sharpness frequency plot
nexttile(7);
x3 = categorical({'Blurry', 'Moderate', 'Sharp'});
x3 = reordercats(x3,{'Blurry', 'Moderate', 'Sharp'});
bc3 = bar(x3, ff2(7:9), 'FaceColor', 'flat');
bc3.CData = [cm3(10,:); cm3(110,:); cm3(210,:)];
hold on
eb3 = errorbar(x3, ff2(7:9), err(7:9), '.');
eb3.Color = [0.1 0.1 0.1];
eb3.CapSize = 5;
eb3.MarkerSize = 0.1;
box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 12,...
    'TickLength', [0.03 0.03])
ylim([0 min(1.1 * (max(ff2(7:9) + err(7:9))), 100)])
xlabel('Sharpness level', 'FontName', 'SansSerif', 'FontSize', 14,...
    'FontWeight', 'bold')
ylabel('Frequency (%)', 'FontName', 'SansSerif', 'FontSize', 14,...
    'FontWeight', 'bold')
hold off

% morphological type frequency plot
nexttile(2, [1,2]);
x4 = categorical({'Fractal soot', 'Dense soot', 'Coated soot',...
    'Tarball', 'Greyball', 'Hybrid', 'Miscellaneous'});
x4 = reordercats(x4,{'Fractal soot', 'Dense soot', 'Coated soot',...
    'Tarball', 'Greyball', 'Hybrid', 'Miscellaneous'});
bc4 = bar(x4, ff2(10:16), 'FaceColor', 'flat');
bc4.CData = [cm4(20,:); cm4(50,:); cm4(80,:); cm4(110,:);...
    cm4(140,:); cm4(170,:); cm4(200,:)];
hold on
eb4 = errorbar(x4, ff2(10:16), err(10:16), '.');
eb4.Color = [0.1 0.1 0.1];
eb4.CapSize = 5;
eb4.MarkerSize = 0.1;
box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 12,...
    'TickLength', [0.01 0.01])
ylim([0 min(1.1 * (max(ff2(10:16) + err(10:16))), 100)])
xlabel('Morphological type', 'FontName', 'SansSerif', 'FontSize', 14,...
    'FontWeight', 'bold')
ylabel('Frequency (%)', 'FontName', 'SansSerif', 'FontSize', 14,...
    'FontWeight', 'bold')

title(tt, 'Visually interpreted morphological properties',...
    'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 16)
hold off

if nargout < 4
    clear hf % delete the results figure handle if not requested as an...
        % ...output
end

% size frequency plot
nexttile(5, [2,2]);
da = cat(1, Aggs.da); % area equivalent diameter
x5_edg = min(da) + (max(da) - min(da)) * log(1 : (exp(1) - 1)/10 : exp(1));
    % lograthmic size binning
his_da = histogram(da, x5_edg); % frequency histogram of sizes
x5 = zeros(10,1); % bin middle points
for i = 1 : 10
    x5(i) = sqrt(his_da.BinEdges(i) * his_da.BinEdges(i+1));
end
y5 = his_da.Values;
bc5 = bar(x5, 100 * y5 / n_agg_tot, 'BarWidth', 1, 'FaceColor', 'flat');
i_bc5 = round(1 + (length(cm5) - 1) .* (0 : 1 / 9 : 1));
bc5.CData = cm5(i_bc5,:);
hold on
err_da = 100 * sqrt(y5 / n_agg_tot^2 + y5.^2 / n_agg_tot^3) ./ y5;
    % error propagation of size bins
eb5 = errorbar(x5, 100 * y5 / n_agg_tot, err_da, '.');
eb5.Color = [0.1 0.1 0.1];
eb5.CapSize = 5;
eb5.MarkerSize = 0.1;
box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 12,...
    'TickLength', [0.01 0.01])
set(gca, 'XScale', 'log')
ylim([0 min(1.1 * (max(100 * y5 / n_agg_tot + err_da)), 100)])
xlabel('Projected area-equivalent diameter (nm)', 'FontName',...
    'SansSerif', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('Frequency (%)', 'FontName', 'SansSerif', 'FontSize', 14,...
    'FontWeight', 'bold')
hold off

title(tt, 'User interpreted morphological properties',...
    'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 16)

if nargout < 4
    clear hf % delete the results figure handle if not requested as an...
        % ...output
end

end

