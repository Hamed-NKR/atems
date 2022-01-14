function [fn, Aggs, imgs_binary, hf] =...
    categorize_semiauto(imgs, n_bin, pixsizes, opts)

% ANALYZE_MORPH characterizes the aggregates based on their morphological...
%   ...properties including circularity and optical depth
% ----------------------------------------------------------------------- %
% imgs: cell array of "cropped" grayscale images to be analyzed
% n_bin: number of post-processing bins (4 member vector for size,...
%   ...circularity, optical depth, and sharpness; also can be a single...
%   ...number for all)
% pixsizes: image pixel sizes
% opts: options structure
% fn: A structure containing the frequencies data
% Aggs: aggregate properties tables
% imgs_binary: aggregate segmented binary images
% hf: A structure for figure handles of frequency plots
% ----------------------------------------------------------------------- %

%% Part1: Categorize the particles %%

n_imgs = length(imgs); % number of images

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
    opts_sizing = 'manual'; % default to be manual aggregate sizing
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
Aggs = agg.analyze_binary(imgs_binary, pixsizes, imgs);

agg_id = cat(1, Aggs.id); % aggregate ids
img_id = cat(1, Aggs.img_id); % image ids corresponding to the aggregates
n_agg_i = zeros(n_imgs,1); % initialize number of aggregates within each...
    % ...image
n_agg_tot0 = length(Aggs); % total number of aggregates across all images

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

fn0 = zeros(n_agg_tot0,5); % initialize array to store bin labels of...
    % ...different properties

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
                    '\fontsize{12}\rm\it(FS/Fs/fs: Fractal soot, Compact/Coated soot,',...
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
                    fn0(jj,5) = 1;
                    Aggs(jj).Type = 'Fractal soot';
                case {'CS', 'Cs', 'cs'}
                    fn0(jj,5) = 2;
                    Aggs(jj).Type = 'Compact/coated soot';
                case {'TB', 'Tb', 'tb'}
                    fn0(jj,5) = 3;
                    Aggs(jj).Type = 'Tarball';
                case {'SB', 'Sb', 'sb'}
                    fn0(jj,5) = 4;
                    Aggs(jj).Type = 'Softball';
                case {'H', 'h'}
                    fn0(jj,5) = 5;
                    Aggs(jj).Type = 'Hybrid';
                case {'M', 'm'}
                    fn0(jj,5) = 6;
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

%% Part 2: Bin the properties %%

% ask what to do next
quest = questdlg(['Manual segmentation and categorization completed!',...
    newline, 'Would you like to proceed to automatic categorization?'],...
    'Next step', 'Continue', 'Finish', 'Continue');
% handle response
if strcmp(quest, 'Finish')
    % save important data first
    fout = pwd; % default saving address to be in the working directory
    fout = strcat(fout, '\morphout'); % set the output file folder
    if ~exist('fout', 'dir')
        mkdir(fout)
    end
    fout = strcat(fout, '\morphdata.mat'); % set the output file name
    save(fout, 'Aggs', 'imgs_binary'); % save the results in case an...
        % ...error occurs while plotting
    
    % null frequency & plot outputs
    fn = struct();
    hf = struct();
    
    return
end

% compile the properties from the aggregate population
da = cat(1, Aggs.da); % area equivalent diameter
ca = cat(1, Aggs.ca); % area equivalent circularity
od = cat(1, Aggs.zbar_opt); % optical depth
os = rand(n_agg_tot0,1); % cat(1, Aggs.sbar_opt); % optical sharpness

% assign the number of bins if not given
if ~exist('n_bin', 'var') || isempty(n_bin)
    n_bin = [10, 10, 10, 10];
elseif length(n_bin) == 1
    n_bin = repmat(n_bin, 1, 4);
elseif length(n_bin(:)) ~= 4
    error(['Invalid bin nuumber set!', newline,...
        'Should be a 4 member vector'])
end

% range of properties
del_da = [min(da), max(da)];
del_ca = [min(ca), max(ca)];
del_od = [min(od), max(od)];
del_os = [min(os), max(os)];

% discretize the domains of post-processing variables (binning)
da_bin = del_da(1) + (del_da(2) - del_da(1)) .*...
    log(1 : (exp(1) - 1) / n_bin(1) : exp(1));
ca_bin = del_ca(1) + (del_ca(2) - del_ca(1)) .* (0 : 1 / n_bin(2) : 1);
od_bin = del_od(1) + (del_od(2) - del_od(1)) .* (0 : 1 / n_bin(3) : 1);
os_bin = del_os(1) + (del_os(2) - del_os(1)) .* (0 : 1 / n_bin(4) : 1);

% initialize the bin centers
da_c = zeros(n_bin(1),1);
ca_c = zeros(n_bin(2),1);
od_c = zeros(n_bin(3),1);
os_c = zeros(n_bin(4),1);

% Removing diverged data
ii0 = da > 1e4;
jj0 = (ca > 2) | (ca < -1);
kk0 = (od > 2) | (od < -1);
ll0 = (os > 2) | (os < -1);
rmv = ii0 | jj0 | kk0 | ll0;
da(rmv) = [];
ca(rmv) = [];
od(rmv) = [];
os(rmv) = [];
fn0(rmv,:) = [];

n_agg_tot = length(da); % nubmer of aggregates after removing the...
    % ...diverged data

fn1 = {zeros(n_bin(1),1), zeros(n_bin(2),1), zeros(n_bin(3),1),...
    zeros(n_bin(4),1), zeros(6,1)}; % initialize the 1d number frequency...
        % ...array
fn2 = {zeros(n_bin(1), n_bin(2)), zeros(n_bin(1), n_bin(3)),...
    zeros(n_bin(1), n_bin(4)), zeros(n_bin(1), 6),...
    zeros(n_bin(2), n_bin(3)), zeros(n_bin(2), n_bin(4)),...
    zeros(n_bin(2), 6), zeros(n_bin(3), n_bin(4)),...
    zeros(n_bin(3), 6), zeros(n_bin(4), 6),}; % 2d frequencies



%% 1d frequecies %%

% assign the aggregates to the bins and count the frquencies

% size binning
for i = 1 : n_bin(1)
    if i == 1 % The first bin (only upper limit)
        ii = da < da_bin(i+1);
    elseif i < n_bin(1)
        ii = (da >= da_bin(i)) & (da < da_bin(i+1));
    else % the last bin (all the remaining)
        ii = fn0(:,1) == 0;
    end

    da_c(i) = sqrt(da_bin(i) * da_bin(i+1));
    fn0(ii,1) = i; % column 1 to be size bin labels
    fn1{1}(i) = nnz(ii) / n_agg_tot; % binned size freqs. (%)
end

% circularity binning
for j = 1 : n_bin(2)
    if j == 1
        jj = ca < ca_bin(j+1);
    elseif j < n_bin(2)
        jj = (ca >= ca_bin(j)) & (ca < ca_bin(j+1));
    else
        jj = fn0(:,2) == 0;
    end
    
    ca_c(j) = (ca_bin(j) + ca_bin(j+1)) / 2;
    fn0(jj,2) = j; % column 2 to be circularity bin labels    
    fn1{2}(j) = nnz(jj) / n_agg_tot; % circularity freqs. (%)
end

% optical depth binning
for k = 1 : n_bin(3)
    if k == 1
        kk = od < od_bin(k+1);
    elseif k < n_bin(3)
        kk = (od >= od_bin(k)) & (od < od_bin(k+1));
    else
        kk = fn0(:,3) == 0;
    end
    
    od_c(k) = (od_bin(k) + od_bin(k+1)) / 2;
    fn0(kk,3) = k; % column 3 to be optical depth bin labels
    fn1{3}(k) = nnz(kk) / n_agg_tot; % optical depth freqs. (%)
end

% optical sharpness binning
for l = 1 : n_bin(4)
    if l == 1
        ll = os < os_bin(l+1);
    elseif l < n_bin(4)
        ll = (os >= os_bin(l)) & (os < os_bin(l+1));
    else
        ll = fn0(:,4) == 0;
    end
    
    os_c(l) = (os_bin(l) + os_bin(l+1)) / 2;
    fn0(ll,4) = l; % column 4 to be optical sharpness bin labels
    fn1{4}(l) = nnz(ll) / n_agg_tot; % optical sharpness freqs. (%)
end

% morphological type binning
for m = 1 : 6
    fn1{5}(m) = nnz(fn0(:,5) == m) / n_agg_tot;
end

% initialize the 1d frequency plots
figure;
hf1 = gcf;
hf1.Position = [0, 0, 1500, 1000]; % Position and size
set(hf1, 'color', 'white'); % Background color

% set the figure layout
tt1 = tiledlayout(2,3);
tt1.TileSpacing = 'loose';
tt1.Padding = 'loose';

% circularity subplot
nexttile(1);

bc1 = bar(ca_c, 100 * fn1{2}, 'FaceColor', 'flat');
cm11 = colormap(autumn);
jjj = round(1 + (length(cm11) - 1) .* (0.1 : 0.8 / (n_bin(2) - 1) : 0.9)');
cm11 = cm11(jjj,:); % get the descretized colormap
cm11 = cm11(end:1);
bc1.CData = cm11;
hold on

err1 = sqrt((fn1{2} + fn1{2}.^2) / n_agg_tot); % standard error propagation
eb1 = errorbar(ca_c, 100 * fn1{2}, 100 * err1, '.');
eb1.Color = [0.1 0.1 0.1];
eb1.CapSize = 5;
eb1.MarkerSize = 0.1;

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.02 0.02])

xticks(ca_c)
ylim([0 10 * ceil(10 * max(fn1{2} + err1))])
xlabel('Circularity (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('Frequency (%)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
hold off

% optical depth subplot
nexttile(2);

bc2 = bar(od_c, 100 * fn1{3}, 'FaceColor', 'flat');
cm12 = colormap(gray);
kkk = round(1 + (length(cm12) - 1) .* (0.1 : 0.8 / (n_bin(3) - 1) : 0.9)');
cm12 = cm12(kkk,:);
cm12 = cm12(end:1);
bc2.CData = cm12;
hold on

err2 = sqrt((fn1{3} + fn1{3}.^2) / n_agg_tot);
eb2 = errorbar(od_c, 100 * fn1{3}, 100 * err2, '.');
eb2.Color = [0.6350 0.0780 0.1840];
eb2.CapSize = 5;
eb2.MarkerSize = 0.1;

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.02 0.02])
xticks(od_c)
ylim([0 10 * ceil(10 * max(fn1{3} + err2))])
xlabel('Optical depth (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('Frequency (%)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
hold off

% sharpness subplot
nexttile(3);
bc3 = bar(os_c, 100 * fn1{4}, 'FaceColor', 'flat');
cm13 = colormap(winter);
lll = round(1 + (length(cm13) - 1) .* (0.1 : 0.8 / (n_bin(4) - 1) : 0.9)');
cm13 = cm13(lll,:);
bc3.CData = cm13;
hold on

err3 = sqrt((fn1{4} + fn1{4}.^2) / n_agg_tot);
eb3 = errorbar(os_c, 100 * fn1{4}, 100 * err3, '.');
eb3.Color = [0.1 0.1 0.1];
eb3.CapSize = 5;
eb3.MarkerSize = 0.1;

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.02 0.02])
xticks(os_c)
ylim([0 10 * ceil(10 * max(fn1{4} + err3))])
xlabel('Optical sharpness (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('Frequency (%)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
hold off

% morphological type subplot
nexttile(4);
x4 = categorical({'Fractal soot', 'Compact/coated soot',...
    'Tarball', 'Softball', 'Hybrid', 'Miscellaneous'});
x4 = reordercats(x4,{'Fractal soot', 'Compact/coated soot',...
    'Tarball', 'Softball', 'Hybrid', 'Miscellaneous'});
bc4 = bar(x4, 100 * fn1{5}, 'FaceColor', 'flat');
cm14 = colormap(summer);
mmm = round(1 + (length(cm14) - 1) .* (0.1 : 0.8 / 5 : 0.9)');
cm14 = cm14(mmm,:);
bc4.CData = cm14;
hold on

err4 = sqrt((fn1{5} + fn1{5}.^2) / n_agg_tot);
eb4 = errorbar(x4, 100 * fn1{5}, 100 * err4, '.');
eb4.Color = [0.1 0.1 0.1];
eb4.CapSize = 5;
eb4.MarkerSize = 0.1;

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.02 0.02])
ylim([0 10 * ceil(10 * max(fn1{4} + err3))])
xlabel('Morphological type', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('Frequency (%)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')

title(tt1, 'Distributions of different morphological properties',...
    'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 16)
hold off

% size distribution subplot
nexttile(5, [1,2]);

fn11 = zeros(n_bin(1),1); % initialize dn/dlog(da)
err5 = zeros(n_bin(1),1); % initialize the errors
cm15 = colormap(turbo);
iii = round(1 + (length(cm15) - 1) .* (0.1 : 0.8 / (n_bin(1) - 1) : 0.9)');
cm15 = cm15(iii,:);
for i = 1 : n_bin(1)
    fn11(i) = fn1{1}(i) / log(da_bin(i+1) / da_bin(i));
    err5(i) = sqrt((fn1{1}(i) + fn1{1}(i).^2) / n_agg_tot) /...
        log(da_bin(i+1) / da_bin(i));
    
    rectangle('Position',[da_bin(i), 0, (da_bin(i+1) - da_bin(i)),...
        fn11(i)], 'FaceColor', cm15(i,:));
    hold on
end
hold on

eb5 = errorbar(da_c, fn11, err5, '.');
eb5.Color = [0.1 0.1 0.1];
eb5.CapSize = 5;
eb5.MarkerSize = 0.1;

box on
axis padded
set(gca, 'FontName', 'SansSerif', 'FontSize', 11, 'TickLength',...
    [0.02 0.02])
set(gca, 'XScale', 'log')
set(gca, 'Layer', 'top')
xticks(da_bin)
xlim([del_da(1), del_da(2)])
ylim([0, ceil(10 * max(fn11 + err5)) / 10])
xlabel('Projected area-equivalent diameter (nm)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Frequency size distribution (nm^-^1)',...
    'FontName', 'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
hold off

title(tt1, 'Frequency of morphological properties',...
    'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 16)

%% 2d frequencies %%

% binning & counting

% size distribution vs. circularity
for i = 1 : n_bin(1)
    for j = 1 : n_bin(2)
        fn2{1}(i,j) = nnz((fn0(:,1) == i) & (fn0(:,2) == j)) / n_agg_tot;
    end
end

% size distribution vs. optical depth
for i = 1 : n_bin(1)
    for k = 1 : n_bin(3)
        fn2{2}(i,k) = nnz((fn0(:,1) == i) & (fn0(:,3) == k)) / n_agg_tot;
    end
end

% size distribution vs. sharpness
for i = 1 : n_bin(1)
    for l = 1 : n_bin(4)
        fn2{3}(i,l) = nnz((fn0(:,1) == i) & (fn0(:,4) == l)) / n_agg_tot;
    end
end

% size distribution vs. morphological type
for i = 1 : n_bin(1)
    for m = 1 : 6
        fn2{4}(i,m) = nnz((fn0(:,1) == i) & (fn0(:,5) == m)) / n_agg_tot;
    end
end

% circularity vs. optical depth
for j = 1 : n_bin(2)
    for k = 1 : n_bin(3)
        fn2{5}(j,k) = nnz((fn0(:,2) == j) & (fn0(:,3) == k)) / n_agg_tot;
    end
end

% circularity vs. sharpness
for j = 1 : n_bin(2)
    for l = 1 : n_bin(4)
        fn2{6}(j,l) = nnz((fn0(:,2) == j) & (fn0(:,4) == l)) / n_agg_tot;
    end
end

% circularity vs. morphological type
for j = 1 : n_bin(2)
    for m = 1 : 6
        fn2{7}(j,m) = nnz((fn0(:,2) == j) & (fn0(:,5) == m)) / n_agg_tot;
    end
end

% optical depth vs. sharpness
for k = 1 : n_bin(3)
    for l = 1 : n_bin(4)
        fn2{8}(k,l) = nnz((fn0(:,3) == k) & (fn0(:,4) == l)) / n_agg_tot;
    end
end

% optical depth vs. morphological type
for k = 1 : n_bin(3)
    for m = 1 : 6
        fn2{9}(k,m) = nnz((fn0(:,3) == k) & (fn0(:,5) == m)) / n_agg_tot;
    end
end

% sharpness vs. morphological type
for l = 1 : n_bin(4)
    for m = 1 : 6
        fn2{10}(l,m) = nnz((fn0(:,4) == l) & (fn0(:,5) == m)) / n_agg_tot;
    end
end

% plotting the results
figure;
hf2 = gcf;
hf2.Position = [0, 0, 1000, 2000]; % Position and size
set(hf2, 'color', 'white'); % Background color

% set the figure layout
tt2 = tiledlayout(6,6);
tt2.TileSpacing = 'loose';
tt2.Padding = 'loose';

% size distribution vs circularity subplot
nexttile(1, [1,2]);

pc1 = pcolor(da_bin, ca_bin, [(fn2{1})', zeros(n_bin(2), 1);...
    zeros(1, n_bin(1) + 1)]);
colormap(hot)
% pc1.FaceColor = 'interp';
% pc1.EdgeColor = 'none';

box on
axis padded
set(gca, 'XScale', 'log')
xlim(del_da)
ylim(del_ca)
% xticks(da_bin)
% xtickformat('%.0f')
% xtickangle(90)
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.03 0.03])
xlabel('Projected area-equivalent diameter (nm)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Circularity (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
cb1 = colorbar;
cb1.FontName = 'SansSerif';
cb1.FontSize = 10;
hold off

% size distribution vs optical depth subplot
nexttile(13, [1,2]);

pc2 = pcolor(da_bin, od_bin, [(fn2{2})', zeros(n_bin(3), 1);...
    zeros(1, n_bin(1) + 1)]);
colormap(gray)
% pc2.FaceColor = 'interp';
% pc2.EdgeColor = 'none';

box on
axis padded
set(gca, 'XScale', 'log')
xlim(del_da)
ylim(del_od)
% xticks(da_bin)
% xtickformat('%.0f')
% xtickangle(90)
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.03 0.03])
xlabel('Projected area-equivalent diameter (nm)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Optical depth (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
cb2 = colorbar;
cb2.FontName = 'SansSerif';
cb2.FontSize = 10;
hold off

% size distribution vs sharpness subplot
nexttile(25, [1,2]);

pc3 = pcolor(da_bin, os_bin, [(fn2{3})', zeros(n_bin(4), 1);...
    zeros(1, n_bin(1) + 1)]);
colormap(parula)
% pc3.FaceColor = 'interp';
% pc3.EdgeColor = 'none';

box on
axis padded
set(gca, 'XScale', 'log')
xlim(del_da)
ylim(del_os)
% xticks(da_bin)
% xtickformat('%.0f')
% xtickangle(90)
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.03 0.03])
xlabel('Projected area-equivalent diameter (nm)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Optical sharpness (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
cb3 = colorbar;
cb3.FontName = 'SansSerif';
cb3.FontSize = 10;
hold off

% circularity vs optical depth subplot
nexttile(2, [1,2]);

isc4 = imagesc(ca_bin, od_bin, (fn2{5})');
colormap(copper)
% isc4.Interpolation = 'bilinear';

box on
axis padded
xlim(del_ca)
ylim(del_od)
xticks(ca_bin)
xtickformat('%.0f')
yticks(od_bin)
ytickformat('%.0f')
grid on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.03 0.03], 'GridColor', [0, 0, 0], 'GridAlpha', 1)
xlabel('Circularity (-)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Optical depth (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
cb4 = colorbar;
cb4.FontName = 'SansSerif';
cb4.FontSize = 10;
hold off

% circularity vs sharpness subplot
nexttile(14, [1,2]);

isc5 = imagesc(ca_bin, os_bin, (fn2{6})');
colormap(spring)
% isc5.Interpolation = 'bilinear';

box on
axis padded
xlim(del_ca)
ylim(del_os)
xticks(ca_bin)
xtickformat('%.0f')
yticks(os_bin)
ytickformat('%.0f')
grid on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.03 0.03], 'GridColor', [0, 0, 0], 'GridAlpha', 1)
xlabel('Circularity (-)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Optical sharpness (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
cb5 = colorbar;
cb5.FontName = 'SansSerif';
cb5.FontSize = 10;
hold off

% optical depth vs sharpness subplot
nexttile(26, [1,2]);

isc6 = imagesc(od_bin, os_bin, (fn2{8})');
colormap(bone)
% isc6.Interpolation = 'bilinear';

box on
axis padded
xlim(del_od)
ylim(del_os)
xticks(od_bin)
xtickformat('%.0f')
yticks(os_bin)
ytickformat('%.0f')
grid on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.03 0.03], 'GridColor', [0, 0, 0], 'GridAlpha', 1)
xlabel('Optical depth (-)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Optical sharpness (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
cb6 = colorbar;
cb6.FontName = 'SansSerif';
cb6.FontSize = 10;
hold off

% size distribution vs morphlogical type
nexttile(3, [2,3]);

x7 = categorical({'Fractal soot', 'Compact/coated soot',...
    'Tarball', 'Softball', 'Hybrid', 'Miscellaneous', ''});
x7 = reordercats(x7,{'Fractal soot', 'Compact/coated soot',...
    'Tarball', 'Softball', 'Hybrid', 'Miscellaneous', ''});
pc7 = pcolor(x7, da_bin, [(fn2{4})', zeros(6, 1);...
    zeros(1, n_bin(1) + 1)]);
colormap(summer)
% pc7.FaceColor = 'interp';
% pc7.EdgeColor = 'none';

box on
axis padded
set(gca, 'YScale', 'log')
xtickangle(45)
ylim(del_da)
% yticks(da_bin)
% ytickformat('%.0f')
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.03 0.03])
xlabel('Morphological type', 'FontName', 'FontName', 'SansSerif',...
    'FontSize', 12, 'FontWeight', 'bold')
ylabel('Projected area-equivalent diameter (nm)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
cb7 = colorbar;
cb7.FontName = 'SansSerif';
cb7.FontSize = 10;
hold off

% circularity vs morphlogical type
nexttile(21, [2,3]);

isc8 = imagesc(x7, ca_bin, (fn2{7})');
colormap(autumn)
% isc8.Interpolation = 'bilinear';

box on
axis padded
xtickangle(45)
ylim(del_ca)
yticks(ca_bin)
ytickformat('%.0f')
grid on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.03 0.03])
xlabel('Morphological type', 'FontName', 'FontName', 'SansSerif',...
    'FontSize', 12, 'FontWeight', 'bold')
ylabel('Circularity (-)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
cb8 = colorbar;
cb8.FontName = 'SansSerif';
cb8.FontSize = 10;
hold off

% optical depth vs morphlogical type
nexttile(5, [2,3]);

isc9 = imagesc(x7, od_bin, (fn2{9})');
colormap(gray)
% isc9.Interpolation = 'bilinear';

box on
axis padded
xtickangle(45)
ylim(del_od)
yticks(od_bin)
ytickformat('%.0f')
grid on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.03 0.03])
xlabel('Morphological type', 'FontName', 'FontName', 'SansSerif',...
    'FontSize', 12, 'FontWeight', 'bold')
ylabel('Optical depth (-)', 'FontName', 'SansSerif',...
    'FontSize', 12, 'FontWeight', 'bold')
cb9 = colorbar;
cb9.FontName = 'SansSerif';
cb9.FontSize = 10;
hold off

% optical sharpness vs morphlogical type
nexttile(23, [2,3]);

isc10 = imagesc(x7, os_bin, (fn2{10})');
colormap(winter)
% isc10.Interpolation = 'bilinear';

box on
axis padded
xtickangle(45)
ylim(del_os)
yticks(os_bin)
ytickformat('%.0f')
grid on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.03 0.03])
xlabel('Morphological type', 'FontName', 'FontName', 'SansSerif',...
    'FontSize', 12, 'FontWeight', 'bold')
ylabel('Optical sharpness (-)', 'FontName', 'SansSerif',...
    'FontSize', 12, 'FontWeight', 'bold')
cb10 = colorbar;
cb10.FontName = 'SansSerif';
cb10.FontSize = 10;
hold off

title(tt2, 'Bivariate frequency of morphological properties',...
    'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 16)

% parametrical colorcoded scatter plots
figure;
hf3 = gcf;
hf3.Position = [0, 0, 700, 2100]; % Position and size
set(hf3, 'color', 'white'); % Background color

% set the figure layout
tt3 = tiledlayout(2,3);
tt3.TileSpacing = 'loose';
tt3.Padding = 'loose';

% size vs circularity
cm31 = colormap(hot);
m4 = round(1 + (length(cm31) - 1) .* (0.1 : 0.8 / 5 : 0.9)');
cm31 = cm31(m4,:);
m5 = cell(5,1);

nexttile
for m = 1 : 6
    m5{i} = (fn0(:,5) == 5);
    scatter(da(m5{i}), ca(m5{i}), 10, cm31(m,:))
    hold on
end

box on
axis padded
xlim(del_da)
ylim(del_ca)
set(gca, 'XScale', 'log')
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
xlabel('Projected area equivalent diamter (nm)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Circularity (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
hold off

% size vs optical depth
cm32 = colormap(gray);
cm32 = cm32(m4,:);

nexttile
for m = 1 : 6
    scatter(da(m5{i}), od(m5{i}), 10, cm32(m,:))
    hold on
end

box on
axis padded
xlim(del_da)
ylim(del_od)
set(gca, 'XScale', 'log')
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
xlabel('Projected area equivalent diamter (nm)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Optical depth (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
hold off

% size vs sharpness
cm33 = colormap(parula);
cm33 = cm33(m4,:);

nexttile
for m = 1 : 6
    scatter(da(m5{i}), os(m5{i}), 10, cm33(m,:))
    hold on
end

box on
axis padded
xlim(del_da)
ylim(del_os)
set(gca, 'XScale', 'log')
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
xlabel('Projected area equivalent diamter (nm)', 'FontName',...
    'SansSerif', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Optical sharpness (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
hold off

% circularity vs optical depth
cm34 = colormap(copper);
cm34 = cm34(m4,:);

nexttile
for m = 1 : 6
    scatter(ca(m5{i}), od(m5{i}), 10, cm34(m,:))
    hold on
end

box on
axis padded
xlim(del_ca)
ylim(del_od)
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
xlabel('Circularity (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('Optical depth (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
hold off

% circularity vs optical sharpness
cm35 = colormap(copper);
cm35 = cm35(m4,:);

nexttile
for m = 1 : 6
    scatter(ca(m5{i}), os(m5{i}), 10, cm35(m,:))
    hold on
end

box on
axis padded
xlim(del_ca)
ylim(del_os)
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
xlabel('Circularity (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('Optical sharpness (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
hold off

% optical depth vs optical sharpness
cm36 = colormap(copper);
cm36 = cm36(m4,:);

nexttile
for m = 1 : 6
    scatter(od(m5{i}), os(m5{i}), 10, cm36(m,:))
    hold on
end

box on
axis padded
xlim(del_od)
ylim(del_os)
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
xlabel('Optial depth (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('Optical sharpness (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
hold off

legtxt3 = {'Fractal soot', 'Compact/coated soot', 'Tarball', 'Softball',...
    'Hybrid', 'Miscellaneous'};
legend(legtxt3, 'Location', 'eastoutside', 'FontName', 'SansSerif',...
    'FontSize', 12);

%% Store outputs %%

fn = struct('1d', fn1, '2d', fn2, 'labels', fn0, 'bins', {da_bin,...
    ca_bin, od_bin, os_bin, x4}); % generate the output frequency structure

% delete the results figure handles if not requested as outputs;...
    % ...otherwise, compile them in a structure for output
if nargout < 4
    clear hf1 hf2 hf3
else
    hf = struct('1d', hf1, '2d', hf2, '3d', hf3);
end

end
