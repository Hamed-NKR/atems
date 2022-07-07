function [I, G_I] = spoptics0(Aggs, id_im, id_par, opts)
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

if ~exist('opts', 'var')
    opts = struct();
end

% Setting defaults for the sampling line

% number of lines to sample on
if ~isfield(opts, 'n_lin')
    opts.n_lin = [];
end
n_lin = opts.n_lin;
if isempty(n_lin) || (n_lin < 6)
    n_lin = 6;
else
    n_lin = round(n_lin);
end

% number of points to on the sample lines
if ~isfield(opts, 'n_p')
    opts.n_p = [];
end
n_p = opts.n_p;
if isempty(n_p)
    n_p = 1000 * ones(n_lin,1);
elseif length(n_p) ~= n_lin
    n_p = repmat(round(n_p(1)), n_lin, 1);
end

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

% set the color set
lc = colormap(turbo); % initialize the colormap of sampling lines
ii = round(1 + (length(lc) - 1) .* (0.05 : 0.9 / (n_lin - 1) : 0.95)');
lc = lc(ii,:); % descretize the colormap

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

% get the nearest and furthest boundary points
parbnd = bwboundaries(img_binary); % get the edge pixel locations
nb = length(parbnd{1});
xb = parbnd{1}(:,1);
yb = parbnd{1}(:,2);
[xb, yb] = poly2cw(xb, yb); % organize edge points
com = Aggs(id2_par).center_mass; % rown and column indices of agg com
R_b = sqrt((xb - com(1)).^2 + (yb - com(2)).^2); % get the edge point distances from COM
ii1_b = find(R_b == max(R_b), 1);
R1_b = [xb(ii1_b), yb(ii1_b)];
ii2_b = find(R_b == min(R_b), 1);
R2_b = [xb(ii2_b), yb(ii2_b)];

R = [(0.2 : 0.8 / (n_lin - 5) : 1) * R1_b; (R1_b + R2_b) / 2; R2_b;...
    1.1 * R2_b; 1.2 * R2_b]; % sampling line radii
x = cell(n_lin,1); % x coordinate indices of sampling line 
y = cell(n_lin,1); % y coordinates

l_im = size(Aggs(id2_par).image); % total size of the image analyzed

[yn, xn] = meshgrid(1 : l_im(2), 1 : l_im(1)); % neighboring cells coordinates...
    % ...for linear interpolation

I = cell(n_lin,1); % initialize pixel intensity on the sampling lines
G_I = cell(n_lin,1); % initialize gradient of pixel intensity
[Gx_I0, Gy_I0] = gradient(Aggs(id2_im).image); % Gradient of the image
G_I0 = sqrt(Gx_I0.^2 + Gy_I0.^2); % magnitude of gradient

for i = 1 : n_lin
    % angle increments for finding coordinates
    theta = 0 : 2 * pi / n_p(i) : 2 * pi;
    theta = theta(1 : end - 1);
    
    % get the line points
    x{i} = com(1) +  R(i) * sin(theta');
    y{i} = com(2) +  R(i) * cos(theta');
    
    % correct for outside image points
    x{i}(x{i} < 1) = 1;
    x{i}(x{i} > l_im(1)) = l_im(1);
    y{i}(y{i} < 1) = 1;
    y{i}(y{i} > l_im(2)) = l_im(2);
    
    % plot the line
    plot(y{i}, x{i}, 'Color', lc(i,:), 'LineWidth', 2)

    
%     % neighboring cells coordinates
%     xn = [floor(x{i}), ceil(x{i}), floor(x{i}), ceil(x{i})];
%     yn = [floor(y{i}), floor(y{i}), ceil(y{i}), ceil(y{i})];
        
    % interpolate the intensity and its gradient
    I{i} = interp2(xn, yn, Aggs(id2_im).image, x{i}, y{i});
    G_I{i} = interp2(xn, yn, G_I0, x{i}, y{i});
end

% plot intensity on the lines
tt3 = nexttile(2);

lgdtex = cell(n_lin, 1); % legend labels placeholder

dr = cell(n_lin, 1); % distance from com in terms of pixels

for i = 1 : n_lin
    dr{i} = sqrt((x{i} - com(1)).^2 + (y{i} - com(2)).^2);
    
    plot(dr{i}, I{i}, 'Color', lc(i,:), 'LineWidth', 2)
    hold on
    
    % make legends
    if i < (n_lin - 4)
        lgdtex{i} = strcat('R / R_{b_i} =', {' '}, num2str(R{i} / R1_b, '%.2f'));
    elseif i == (n_lin - 4)
        lgdtex{i} = 'R = (R_{b_i} + R_{b_o}) / 2';
    else
        lgdtex{i} = strcat('R / R_{b_o} =', {' '}, num2str(R{i} / R2_b, '%.2f'));
    end
end

box on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out')
xlabel('r (pix)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('I (-)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
lgd = legend(cat(2,lgdtex{:}), 'interpreter', 'tex', 'FontName', 'SansSerif',...
    'FontSize', 10, 'Orientation', 'Horizontal');
lgd.Layout.Tile = 'south';
title(tt3, 'Intensity distrubutions')
hold off

tt4 = nexttile(4);
for i = 1 : n_lin
    plot(dr{i}, G_I{i}, 'Color', lc(i,:), 'LineWidth', 2)
    hold on
end
box on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out')
xlabel('r (pix)', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
ylabel('dI / dr (pix^{-1})', 'FontName', 'SansSerif', 'FontSize', 12,...
    'FontWeight', 'bold')
title(tt4, 'Intensity gradient distrubution')
hold off

