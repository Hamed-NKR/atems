function [I, Ic, GI, GIc] = spoptics(Aggs, id_im, id_par, opts)
% SPOPTICS tracks the variations of optical properties within single...
%   ...particles.
% ----------------------------------------------------------------------- %
% Aggs: particle table of properties
% id_im: index of image to be analyze
% id_par: particle index within the image
% opts: a strcuture containing the options available for analysis and plotting
% I: image intensity distribution
% Ic: cumulative intensity distribution
% GI: distribution of intensity gradient magnitude
% GIc: cumulative distribution of intensity gradient magnitude
% ----------------------------------------------------------------------- %

if ~exist('opts', 'var')
    opts = struct();
end

%%% Setting defaults for the sampling line %%%

% domain extension factor for sampling
if ~isfield(opts, 'cc')
    opts.cc = [];
end
cc = opts.cc;
if isempty(cc) || (cc <= 1)
    cc = 1.2;
    % warning('Unsuitable extention factor; automatically adjusted...')
end

% number of radial points used for the distribution sampling
if ~isfield(opts, 'nr')
    opts.nr = [];
end
nr = opts.nr;
if isempty(nr)  || nr < 50
    nr = 100;
else
    nr = round(nr(1));
end

% number of angular increments for distribution averaging
if ~isfield(opts, 'na')
    opts.na = [];
end
na = opts.na;
if isempty(na) || na < 20
    na = 360;
else
    na = round(na(1));
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
h.Position = [0, 0, 1100, 900]; % position and size
set(h, 'color', 'white'); % background color

% set the figure layout
tt = tiledlayout(2,2);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

% set the color set
lc = colormap(hot); % initialize the colormap of sampling lines
ii = round(1 + (length(lc) - 1) .* (0.05 : 0.9 / 7 : 0.95)');
lc = lc(ii,:); % descretize the colormap

% original image
tt1 = nexttile(1);
id2_im = im_ids(id_im);
tools.imshow(Aggs(id2_im).image)
title(tt1, 'Original TEM image', 'FontSize', 14)

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
title(tt2, 'Particle selected & sampling lines', 'FontSize', 14)
hold on

% get the nearest and furthest boundary points
parbnd = bwboundaries(Aggs(id2_par).binary); % get the edge pixel locations
% nb = length(parbnd{1});
xb = parbnd{1}(:,1);
yb = parbnd{1}(:,2);
[xb, yb] = poly2cw(xb, yb); % organize edge points
com = Aggs(id2_par).center_mass; % rown and column indices of agg com
Rb = sqrt((xb - com(1)).^2 + (yb - com(2)).^2); % get the edge point distances from COM
R1b = min(Rb);
R2b = max(Rb);

R0 = [(0.25 : 0.25 : 1)' * R1b; (R1b + R2b) / 2; R2b;...
    1.1 * R2b; 1.2 * R2b]; % radii of highlighted sampling lines for visualization
R = (0 : 1 / (nr - 1) : 1)' * cc * R2b; % descretized radii for distribution sampling 

l_im = size(Aggs(id2_im).image); % overall size of the image analyzed

[yn, xn] = meshgrid(1 : l_im(2), 1 : l_im(1)); % generate an index map for later intensity averaging

% initialize placeholders for the distributions to be sampled
I = zeros(nr,1); 
GI = zeros(nr,1);
Ic = zeros(nr,1); 
GIc = zeros(nr,1);

I0 = double(Aggs(id2_im).image); % image intensity
[GxI0, GyI0] = gradient(I0); % gradients of the image
GI0 = sqrt(GxI0.^2 + GyI0.^2); % magnitude of gradient throughout the image

% angle increments for finding sample line coordinates
theta = 0 : 2 * pi / na : 2 * pi;
theta = theta(1 : end - 1);

% plot the highlight circles
phi = 0 : pi / 180 : 2 * pi;
for ii = 1 : 8
    x0 = com(1) +  R0(ii) * sin(phi');
    y0 = com(2) +  R0(ii) * cos(phi');
    
    plot(y0, x0, 'Color', lc(ii,:), 'LineWidth', 1)
end
hold off

disp(' ')
disp('Pefroming radial intensity sampling...')
tools.textbar([0, nr]); % initialize textbar

for i = 1 : nr    
    % get the line points
    x = com(1) +  R(i) * sin(theta');
    y = com(2) +  R(i) * cos(theta');
    
    % correct for outside image points
    x(x < 1) = 1;
    x(x > l_im(1)) = l_im(1);
    y(y < 1) = 1;
    y(y > l_im(2)) = l_im(2);
    
    % interpolate the intensity and its gradient, take average and accumulate
    I(i) = mean(interp2(xn', yn', I0', x, y));
    Ic(i) = sum(I) / nnz(I);
    GI(i) = mean(interp2(xn', yn', GI0', x, y));
    GIc(i) = sum(GI) / nnz(GI);
    
    tools.textbar([i, nr]); % update textbar
end

% generate legend labels
lgdtxt = {'R = 0.25 R_{b_i}', 'R = 0.5 R_{b_i}', 'R = 0.75 R_{b_i}',...
    'R = R_{b_i}', 'R = (R_{b_i} + R_{b_o}) / 2', 'R = R_{b_o}',...
    'R = 1.1 R_{b_o}', 'R = 1.2 R_{b_o}'}; 
lgd = legend(lgdtxt, 'interpreter', 'tex', 'FontName', 'SansSerif',...
    'FontSize', 11);
lgd.Layout.Tile = 'east';
hold off

% plot radial intensity distribution
tt3 = nexttile(2);

yyaxis left
plot(R, I, 'LineWidth', 2);
ylabel('Local mean intensity (-)', 'FontName', 'SansSerif', 'FontSize', 12)
hold on

yyaxis right
plot(R, Ic, 'LineWidth', 2);
ylabel('Cumulative mean intensity (-)', 'FontName', 'SansSerif', 'FontSize', 12)

y3 = ylim;
for ii = 1 : 8
    line([R0(ii), R0(ii)], [y3(1), y3(2)], 'Color', lc(ii,:),...
        'LineWidth', 1.5, 'LineStyle', '--')
end

box on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out')
xlabel('Radial distance from COM (pix)', 'FontName', 'SansSerif', 'FontSize', 12)
xlim([R(1), R(end)])
title(tt3, 'Intensity distrubution', 'FontSize', 14)
hold off

% plot the gradient of radial intensity distribution
tt4 = nexttile(4);

yyaxis left
plot(R, GI, 'LineWidth', 2);
ylabel('Local mean intensity gradient (pix^{-1})', 'FontName', 'SansSerif', 'FontSize', 12)
hold on

yyaxis right
plot(R, GIc, 'LineWidth', 2);
ylabel('cumulative mean intensity gradient (pix^{-1})', 'FontName', 'SansSerif', 'FontSize', 12)

y4 = ylim;
for ii = 1 : 8
    line([R0(ii), R0(ii)], [y4(1), y4(2)], 'Color', lc(ii,:),...
        'LineWidth', 1.5, 'LineStyle', '--')
end

box on
set(gca, 'FontName', 'SansSerif', 'FontSize', 11,...
    'TickLength', [0.01 0.01], 'TickDir', 'out')
xlabel('r (pix)', 'FontName', 'SansSerif', 'FontSize', 12)
xlim([R(1), R(end)])
title(tt4, 'Intensity gradient distrubution', 'FontSize', 14)
hold off


