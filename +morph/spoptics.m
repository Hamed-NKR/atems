function [I, G_I, h] = spoptics(Aggs, id_im, id_par, opts)
% SPMONITOR tracks the variations of optical properties properties...
%   ...within single particles.
% ----------------------------------------------------------------------- %
% Aggs: aggregates table of properties
% id_im: index of image to be analyze
% id_par: particle index within the image
% opts: a strcuture containing the options available for plotting.
% h: output figure handle
% pardat: sampling line optica
% ----------------------------------------------------------------------- %

% initialize options if not given
if ~(exist('opts', 'var') && isfield(opts, 'n_ang') &&...
        isfield(opts, 'c_dil'))
    opts = struct('n_ang', [], 'c_dil', []);
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
if id_im > n_im
    id_im = n_im;
    warning('Image id exceeding the maximum value; adjusted to maximum...')
elseif id_im < 1
    id_im = 1;
    warning('Image id below the minimum value; adjusted to minimum...')
end
if id_par > par_n(id_im)
    id_par = par_n(id_im);
    warning('Particle id exceeding the maximum value; adjusted to maximum...')
elseif id_par < 1
    id_par = 1;
    warning('Particle id below the minimum value; adjusted to minimum...')
end

% initialize figure properties
figure;
h = gcf;
h.Position = [0, 0, 900, 1200]; % position and size
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
se2 = strel('disk',5);
edg_dilated2 = imdilate(par_edg, se2);
com = round(Aggs(id2_par).center_mass);
lin = cell(n_ang,1);
r = cell(n_ang,1);
angs = 2 * pi * (0 : 1 / n_ang : 1);
angs = angs(1 : end - 1);

[x_e, y_e] = find(edg_dilated2 == 1);
[x_e, y_e] = poly2cw(x_e, y_e);

for i = 1 : n_ang
    ii = max(((y_e - com(2)).^2 ./ (x_e - com(1)).^2) - tan(angs(i)) <...
        pi / 180);
    drawline('Position', [com(1) com(2); x_e(ii) y_e(ii)], 'Color', 'r')
    lin{i} = unique(round([(0 : 100)' * (x_e(ii) - com(1)),...
        (0 : 100)' * (y_e(ii) - com(2))]), 'row');
    r{i} = sum(sqrt(lin{i}),2);
    lin{i} = lin{i} + [com(1), com(2)];
end

% pixel intensity on the lines drawn
tt3 = nexttile(2);
I = cell(n_ang,1);
for i = 1 : n_ang
    I{i} = zeros(length(r{i}),1);
    for j = 1 : length(r{i})
        I{i}(j) = Aggs(id2_im).image(lin{i}(j,1), lin{i}(j,2));
        plot(r{i}(j), I{i}(j))
        hold on
    end
end
box on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out')
set(gca, 'XMinorTick', 'off')
xlabel('r (pix)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('I (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
title(tt3, 'Intensity distrubution')
hold off

% gradient of pixel intensity
tt4 = nexttile(4);
G_I = cell(n_ang,1);
for i = 1 : n_ang
    G_I{i} = gradient(I{i});
    plot(r{i}, G_I{i})
    hold on
end
box on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out')
set(gca, 'XMinorTick', 'off')
xlabel('r (pix)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('I (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
title(tt4, 'Intensity gradient distrubution')
hold off

