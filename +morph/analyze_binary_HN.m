
% ANALYZE_BINARY  Label and analyze a binary mask to quantify aggregates.
%  Properties computed include the radius of gyration, projected
%  area-equivalent diameter, aspect ratio, etc.
% 
%  AGGS = analyze_binary(IMGS_BINARY) uses the binary mask to identify
%  independent aggregates in the image. Applies a pixel size of 1 nm/pixel,
%  such that results will be given pixels rather than in nm. The output is 
%  a data structure with one entry per identified aggregate. 
%  NOT RECOMMENDED, as orig. image information is not transferred to 
%  AGGS for subsequent analysis. 
%  
%  AGGS = analyze_binary(IMGS_BINARY,IMGS) uses the IMGS data structure to
%  extract addition properties, such as the cropped image, filename, pixel
%  size.
% 
%  AGGS = analyze_binary(IMGS_BINARY,PIXSIZE) uses the pixel size of each
%  image to determine physical quantities from the aggregates, such as the
%  projected area-equivalent diameter in nm. 
%  NOT RECOMMENDED. Preferred use includes IMGS. 
% 
%  AGGS = analyze_binary(IMGS_BINARY,PIXSIZE,IMGS) uses the IMGS cell array
%  containing the original images in plotting and stores the result for
%  subsequent analysis. 
% 
%  AGGS = analyze_binary(IMGS_BINARY,PIXSIZE,IMGS,FNAME) also adds the
%  filename, FNAME, to the AGGS structure. 
% 
%  AGGS = analyze_binary(IMGS_BINARY,PIXSIZE,IMGS,FNAME,F_EDGES) adds a 
%  flag to determine whether aggregate at the border are cleared. 
%  F_EDGES = 1 removes the border aggregates. 
% 
%  AGGS = analyze_binary(IMGS_BINARY,PIXSIZE,IMGS,FNAME,F_EDGES,F_PLOT) 
%  adds a flag of whether to plot the results as they are computed.
% 
%  AUTHOR: Timothy Sipkens, 2019-11-26

function [Aggs] = analyze_binary(imgs_binary, pixsize, ...
    imgs, fname, f_edges, f_plot)

%-- Parse inputs ---------------------------------------------------------%
if isstruct(imgs_binary) % consider case that structure is given as input
    Imgs = imgs;
    imgs = {Imgs.cropped};
    pixsize = [Imgs.pixsize];
    fname = {Imgs.fname};
else
    % consider case that a single image is given
    if ~iscell(imgs_binary); imgs_binary = {imgs_binary}; end
end

if ~exist('pixsize','var'); pixsize = []; end
if isempty(pixsize); pixsize = ones(size(imgs_binary)); end

if ~exist('imgs','var'); imgs = []; end
if isempty(imgs)
    imgs = imgs_binary;
    for ii=1:length(imgs)
        imgs{ii} = uint8(155 .* ~imgs{ii} + 100);
    end
end
if ~iscell(imgs); imgs = {imgs}; end

if ~exist('fname','var'); fname = []; end
if isempty(fname); fname = cell(size(imgs_binary)); end

% Flag for whether to remove aggregates at the edges of the images.
% 1 removes border aggregates, other values keep border aggregates.
% Default is to remove border aggregates. 
if ~exist('f_edges','var'); f_edges = []; end
if isempty(f_edges); f_edges = 1; end

% Flag for whether to show progress in a figure.
if ~exist('f_plot','var'); f_plot = []; end
if isempty(f_plot); f_plot = 1; end
%-------------------------------------------------------------------------%


if f_plot==1; f0 = figure; end % intialize a new figure to show progress
Aggs = struct([]); % initialize Aggs structure
id = 0;

tools.textheader('Analyzing binaries');
disp(' Progress:'); tools.textbar([0, length(imgs_binary)]);
for ii=1:length(imgs_binary) % loop through provided images
    
    img_binary = imgs_binary{ii};
    img = imgs{ii};
    
    
    % If more than 25% of the image is boundary aggregate, the method likely failed.
    % Skip this image and continue on. This is done before remove border
    % aggregates, if relevant. 
    bwborder = img_binary - imclearborder(img_binary);
    if (nnz(bwborder) / numel(img_binary)) > 0.25
        continue;
    end
    
    
    % Check if any of the borders are >20% aggregate. 
    % This is likely a problem image, try to mask result. 
    nn = [nnz(img_binary(:, 1)) / size(img_binary, 1), ...
        nnz(img_binary(:, end)) / size(img_binary, 1), ...
        nnz(img_binary(1, :)) / size(img_binary, 2), ...
        nnz(img_binary(end, :)) / size(img_binary, 2)];
    if any(nn > 0.2)
        ia = img_binary;
        % Zero edges that are not problematic
        % (avoids real aggreagates on border).
        if ~(nn(1) > 0.2); ia(:, 1) = 0; end
        if ~(nn(2) > 0.2); ia(:, end) = 0; end
        if ~(nn(3) > 0.2); ia(1, :) = 0; end
        if ~(nn(4) > 0.2); ia(end, :) = 0; end
        
        % Examine just remaining border regions.
        img_edges = ia - imclearborder(ia);
        
        % For generating a cirle excluding border regions.
        [i1, i2] = meshgrid(1:size(img_binary, 2), 1:size(img_binary, 1));
        fun = @(x) sqrt((i1 - x(1)).^2 + (i2 - x(2)).^2) > x(3);
        
        img_edm = bwdist(img_edges);  % euclidean distance map (EDM)
        [~, m1] = max(max(img_edm, [], 1));  % find peak of EDM
        [~, m2] = max(max(img_edm, [], 2));
        img_binary = ...
            or(fun([m1, m2, img_edm(m2, m1) - 100]), img_binary);
    end
    
    
    % If clearing aggregate borders...
    if f_edges
        img_binary = imclearborder(img_binary);
    end
    
    
    % Remove aggregates below 10 pixels, which will
    % cause problems with primary particle sizing.
    % Segmentation technqiues may impose different minimums.
    img_binary = bwareaopen(img_binary, 10);
    
    
    % Detect distinct aggregates.
    CC = bwconncomp(img_binary); % find seperate aggregates
    naggs = CC.NumObjects; % count number of aggregates
    
    % Background optical depth
    CC_PixelIdxList = cell2mat(reshape(CC.PixelIdxList, naggs, 1));
    seg_backgr = img(:);
    seg_backgr(CC_PixelIdxList) = []; % segmented background grayscale image
    I_backgr = mean(seg_backgr); % background image intensity
    
    % If more than 50 aggregates were found, the method likely failed. 
    % Skip this image and continue on. 
    if naggs>50; continue; end
    
    % If no aggregates, skip image. 
    if naggs==0; continue; end
    
    
    Aggs0 = struct([]); % re-initialize Aggs0 structure
    Aggs0(naggs).fname = '';
        % pre-allocate new space for aggregates and assign filename

    %== Main loop to analyze each aggregate ==============================%
    if f_plot==1; set(groot,'CurrentFigure',f0); tools.imshow_binary(img, img_binary); end
    
    for jj = 1:naggs % loop through number of found aggregates
        
        id = id + 1; % increment global index counter
        Aggs0(jj).id = id;
        Aggs0(jj).img_id = ii;
        
        Aggs0(jj).fname = fname{ii};
        Aggs0(jj).pixsize = pixsize(ii);
        
        if jj==1
            Aggs0(1).image = img; % store the overall image for the first aggregate
        end
        
        %-- Step 3-2: Prepare an image of the isolated aggregate ---------%
        img_binary = zeros(size(img_binary));
        img_binary(CC.PixelIdxList{1,jj}) = 1;
        Aggs0(jj).binary = logical(img_binary); % store binary image
        
        % Get a cropped version of the aggregate
        % 'autocrop' method included below.
        [~, ~, Aggs0(jj).rect] = autocrop(img, img_binary);
        
        
        %== Compute aggregate dimensions/parameters ======================%
        SE = strel('disk', 1);
        img_dilated = imdilate(img_binary,SE);
        img_edg = img_dilated - img_binary;
        
        [row, col] = find(imcrop(Aggs0(jj).binary, Aggs0(jj).rect));
        Aggs0(jj).length = max((max(row)-min(row)), (max(col)-min(col))) * pixsize(ii);
        Aggs0(jj).width = min((max(row)-min(row)), (max(col)-min(col))) * pixsize(ii);
        Aggs0(jj).aspect_ratio = Aggs0(jj).length / Aggs0(jj).width;
        
        Aggs0(jj).num_pixels = nnz(img_binary); % number of non-zero pixels
        Aggs0(jj).da = ((Aggs0(jj).num_pixels/pi)^.5) * ...
            2 * pixsize(ii); % area-equialent diameter [nm]
        Aggs0(jj).area = nnz(img_binary) .* ... 
            pixsize(ii) ^ 2; % aggregate area [nm^2]
        Aggs0(jj).Rg = gyration(img_binary, pixsize(ii)); % calculate radius of gyration [nm]
        
%         Aggs0(jj).perimeter0 = pixsize(ii) * (sum(sum(img_edg~=0))); % calculate aggregate perimeter
        
        Aggs0(jj).perimeter = pixsize(ii) * get_perimeter2(img_edg);
        
        [x,y] = find(img_binary ~= 0);
        Aggs0(jj).center_mass = [mean(x); mean(y)];
        
        p_c = 2 * sqrt(pi * Aggs0(jj).area); % perimeter of aggregate area-equivalent circle
        
        Aggs0(jj).ca = p_c / Aggs0(jj).perimeter; % aggregate area-equiv. circulirity
            % the degree of being far from a circle (1: circle, 0: straight line)
        
        agg_grayscale = img(CC.PixelIdxList{1,jj}); % the selected agg's grayscale pixel values
            Aggs0(jj).zbar_opt = (I_backgr - mean(agg_grayscale)) / I_backgr; % agg's optical depth metric (1: black, 0: white)
        
        
        
        if f_plot==1; set(groot,'CurrentFigure',f0); tools.imshow_agg(Aggs0, ii, 0); title(num2str(ii)); drawnow; end
    end

    
    if f_plot==1; pause(0.05); end % pause very briefly to show overall aggregates
    
    Aggs = [Aggs, Aggs0]; % append current aggregate data
    
    tools.textbar([ii, length(imgs_binary)]);
end

if f_plot==1; close(f0); end % close figure showing progress

tools.textheader();

end




%== GYRATION =============================================================%
%   Gyration calculates radius of gyration by assuming every pixel as an area
%   of pixsize^2
%   Author: Ramin Dastanpour
%   
% Output:
%   Rg      Radius of gyration [nm]
function [Rg] = gyration(img_binary,pixsize)


total_area = nnz(img_binary)*pixsize^2;

[xpos,ypos] = find(img_binary);
n_pix = size(xpos,1);
Centroid.x = sum(xpos)/n_pix;
Centroid.y = sum(ypos)/n_pix;

Ar2 = zeros(n_pix,1);

for kk = 1:n_pix
    Ar2(kk,1) = (((xpos(kk,1)-Centroid.x)*pixsize)^2+...
        ((ypos(kk,1)-Centroid.y)*pixsize)^2)*pixsize^2;
end

Rg = (sum(Ar2)/total_area)^0.5;

end




%== AUTOCROP =============================================================%
%   Automatically crops an image based on binary information
%   Author:  Yeshun (Samuel) Ma, Timothy Sipkens, 2019-07-23
function [img_cropped,img_binary,rect] = autocrop(img_orig, img_binary)

[x,y] = find(img_binary);

space = 3;
size_img = size(img_orig);

% Find coordinates of top and bottom of aggregate
x_top = min(max(x)+space,size_img(1)); 
x_bottom = max(min(x)-space,1);
y_top = min(max(y)+space,size_img(2)); 
y_bottom = max(min(y)-space,1);


img_binary = img_binary(x_bottom:x_top,y_bottom:y_top);
img_cropped = img_orig(x_bottom:x_top,y_bottom:y_top);
rect = [y_bottom,x_bottom,(y_top-y_bottom),(x_top-x_bottom)];

end

%== GET_PERIMETER2 =============================================================%
%   An updated method to get the perimeter of the aggregate.
function p = get_perimeter2(img_binary)
 
mb = bwboundaries(img_binary);
n_mb = length(mb{1});
 
x_mb = mb{1}(:,2);
y_mb = mb{1}(:,1);
ii_mb = ones(n_mb,1);

[x_mb, y_mb] = poly2cw(x_mb, y_mb);

ii = 1;
for i = 2 : n_mb
    ii_mb(i) = ii;
    if (x_mb(i) ~= x_mb(i-1)) && (y_mb(i) ~= y_mb(i-1))
        ii_mb(i) = ii_mb(i) + 1;
        ii = ii + 1;
    end
end

if (x_mb(1) == x_mb(end)) || (y_mb(1) == y_mb(end))
    ii_mb(ii_mb == ii_mb(end)) = 1;
end

nn_mb = max(ii_mb);
xx_mb = zeros(nn_mb,1);
yy_mb = zeros(nn_mb,1);
p = 0;

for i = 1 : nn_mb
    xx_mb(i) = mean(x_mb(ii_mb == i));
    yy_mb(i) = mean(y_mb(ii_mb == i));
    
    if i > 1
        p = p + sqrt((xx_mb(i) - xx_mb(i-1))^2 +...
           (yy_mb(i) - yy_mb(i-1))^2);
    end
end
p = p + sqrt((xx_mb(1) - xx_mb(end))^2 +...
    (yy_mb(1) - yy_mb(end))^2);
 
end
