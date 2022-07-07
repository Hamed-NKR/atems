function [I, G_I] = spoptics00(Aggs, id_im, id_par, opts)
% SPOPTICS tracks the variations of optical properties within single...
%   ...particles.
% ----------------------------------------------------------------------- %
% Aggs: particle table of properties
% id_im: index of image to be analyze
% id_par: particle index within the image
% opts: a strcuture containing the options available for plotting.
% I: image intensity on the sampling lines
% G_I: intensity gradients on the sampling lines
% ----------------------------------------------------------------------- %

% Setting defaults for the sampling line
if ~exist('opts', 'var')
    opts = struct();
end

if ~isfield(opts, 'n_ang')
    opts.n_ang = [];
end

if ~isfield(opts, 'c_dil')
    opts.c_dil = [];
end

n_ang = opts.n_ang; % number of angles to be analyzed
c_dil = opts.c_dil; % dilation factor for data sampling line

% set defaults for opts
if isempty(n_ang)
    n_ang = 4;
end
if isempty(c_dil)
    c_dil = 3;
end

% adjust the ids if not acceptable
% n_par = length(Aggs); % total particle number
n_im = max(cat(1, Aggs.img_id)); % total image number
im_ids = zeros(n_im,1); % particle ids corresponding to the start of each image
par_n = zeros(n_im,1); % number of particles within each image
Aggs_im_id = cat(1, Aggs.img_id);
for i = 1 : n_im
    im_ids(i) = find(Aggs_im_id == i, 1);
    par_n(i) = nnz(Aggs_im_id == i);
end

% adjust inputs if not acceptable
if id_im > n_im
    id_im = n_im;
    warning('Image id exceeding the maximum value; adjusted to maximum...')
elseif id_im < 1
    id_im = 1;
    warning('Image id below the minimum value; adjusted to minimum...')
elseif mod(id_im, 1) ~= 0
    id_im = round(id_im);
    warning('Image id needs to be integer; input value rounded...')
end
if id_par > par_n(id_im)
    id_par = par_n(id_im);
    warning('Particle id exceeding the maximum value; adjusted to maximum...')
elseif id_par < 1
    id_par = 1;
    warning('Particle id below the minimum value; adjusted to minimum...')
elseif mod(id_par, 1) ~= 0
    id_par = round(id_par);
    warning('Particle id needs to be integer; input value rounded...')
end

% initialize figure properties
figure;
h = gcf;
h.Position = [0, 0, 800, 800]; % position and size
set(h, 'color', 'white'); % background color

% set the figure layout
tt = tiledlayout(2,2);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

% original image
tt1 = nexttile(1);
id2_im = im_ids(id_im);
tools.imshow(Aggs(id2_im).image)
title(tt1, 'Original TEM image')

% image highlighting the selected aggregate and data sampling lines
tt2 = nexttile(3);
id2_par = id2_im + id_par - 1;
cmap = ones(max(max(Aggs(id2_par).binary)), 1) * [0.12, 0.59, 0.96];
im2 = labeloverlay(Aggs(id2_im).image,...
    Aggs(id2_par).binary, 'Transparency', 0.6, 'Colormap', cmap);
par_edg = edge(Aggs(id2_par).binary, 'sobel');
se = strel('disk',1);
edg_dilated = imdilate(par_edg, se);
im2 = uint8(~edg_dilated) .* im2;
tools.imshow(im2); % display the highlighted aggregate and the edge
title(tt2, 'Particle selected & sampling lines')
hold on

% get sampling line start & end coordinates
se2 = strel('disk', c_dil);
edg_dilated2 = imdilate(par_edg, se2);
com = round(Aggs(id2_par).center_mass); % rown and column indices of agg com
r = cell(n_ang,1); % placeholder for sampling line points
dr = cell(n_ang,1); % line points relative to com
rt = zeros(n_ang,2); % coordinates of tips of the sampling line
angs = 2 * pi * (0 : 1 / n_ang : 1); % angles to be sampled
angs = angs(1 : end - 1); % 2 * pi is redundant

[xe, ye] = find(edg_dilated2 == 1);
[xe, ye] = poly2cw(xe, ye);

for i = 1 : n_ang
    % find the line tip
    ii0 = abs(atan2((ye - com(2)), (xe - com(1))) -...
        (angs(i) - pi)) < (pi / 36);
    dr0 = (xe(ii0) - com(1)).^2 + (ye(ii0) - com(2)).^2;
    ii = find(dr0 == max(dr0), 1);
    rt(i,:) = [xe(ii), ye(ii)];
    
    % plot the line
    drawline('Position', [com(2) com(1); rt(i,2) rt(i,1)], 'Color', 'r')
    
    % store the position data
    r{i} = unique(round([(0 : 0.01 : 1)' * (xe(ii) - com(1)),...
        (0 : 0.01 : 1)' * (ye(ii) - com(2))]), 'row');
    dr{i} = sum(sqrt(r{i}),2) * Aggs(id2_par).pixsize;
    r{i} = r{i} + [com(1), com(2)];
end

tt3 = nexttile(2);

I = cell(n_ang,1); % initialize pixel intensity to be sampled

lgdtex = cell(n_ang, 1); % legend labels placeholder
lc = colormap(tt3, jet); % line colormap
ii = round(1 + (length(lc) - 1) .* (0.05 : 0.9 / (n_ang - 1) : 0.95)');
lc = lc(ii,:); % descretize the colormap

for i = 1 : n_ang
    I{i} = zeros(length(dr{i}),1);
    for j = 1 : length(dr{i})
        I{i}(j) = Aggs(id2_im).image(r{i}(j,1), r{i}(j,2)); % store intensity data
    end
    
    plot(dr{i}, I{i}, 'Color', lc(i,:), 'LineWidth', 2)
    hold on
    
    lgdtex{i} = strcat('theta =', {' '}, num2str(angs(i), '%.2f')); % make legends
end
box on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out')
xlabel('r (nm)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('I (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
lgd = legend(cat(2,lgdtex{:}), 'FontName', 'SansSerif', 'FontSize', 10,...
    'Orientation', 'Horizontal');
lgd.Layout.Tile = 'south';
title(tt3, 'Intensity distrubution')
hold off

tt4 = nexttile(4);
G_I = cell(n_ang,1); % initialize gradient of pixel intensity
for i = 1 : n_ang
    G_I{i} = gradient(I{i});
    plot(dr{i}, G_I{i}, 'Color', lc(i,:), 'LineWidth', 2)
    hold on
end
box on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out')
xlabel('r (nm)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('dI / dr (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
title(tt4, 'Intensity gradient distrubution')
hold off

