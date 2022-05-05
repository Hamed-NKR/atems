clc
clear
close all

n_img = [10000, 10000];
% pixsz = 1.9;

img = ones(n_img(1), n_img(2));

cc = [randperm(n_img(1), 1), randperm(n_img(2), 1)];
rc = min(min(n_img(1) - cc(1), cc(1) - 1),...
    min(n_img(2) - cc(2), cc(2) - 1)) / 3;
ac_math = pi * rc^2;
pc_math = 2 * pi * rc;
fprintf('pc_math = %f\n', pc_math);

for i = 1 : n_img(1)
    for j = 1 : n_img(2)
        if sqrt((i-cc(1))^2 + (j-cc(2))^2) <= rc
            img(i,j) = 0;
        end
    end
end

ac_img = numel(img) - nnz(img);
% ac_img = bwarea(~img);
pc_img = 2 * sqrt(pi * ac_img);
% imshow(img);
fprintf('pc_img = %f\n', pc_img);

figure(1)
SE = strel('disk', 1);
% SE = strel('line', 2, 0);
img_dilated = imdilate(img,SE);
edg_dil = img_dilated - img;
pc_edgdil = nnz(edg_dil);
fprintf('pc_edgdil = %f\n', pc_edgdil);
imshow(edg_dil);
% imov1 = imoverlay(edg_dil, img, 'r');
% imshow(imov1);
hold on

mb = bwboundaries(~img);
n_mb = length(mb{1});
pc_mb1 = n_mb;
fprintf('pc_mb1 = %f\n', pc_mb1);
plot(mb{1}(:,2), mb{1}(:,1),'Color','g')

pc_mb2 = 0;
for i = 1 : n_mb
    if i ~= n_mb
        pc_mb2 = pc_mb2 + sqrt(sum((mb{1}(i+1,:) - mb{1}(i,:)).^2));
    else
        pc_mb2 = pc_mb2 + sqrt(sum((mb{1}(1,:) - mb{1}(i,:)).^2));
    end
end
fprintf('pc_mb2 = %f\n', pc_mb2);

x_mb = mb{1}(:,2);
y_mb = mb{1}(:,1);
ii_mb = ones(n_mb,1);

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
pc_mb3 = 0;

for i = 1 : nn_mb
    xx_mb(i) = mean(x_mb(ii_mb == i));
    yy_mb(i) = mean(y_mb(ii_mb == i));
    
    if i > 1
        pc_mb3 = pc_mb3 + sqrt((xx_mb(i) - xx_mb(i-1))^2 +...
           (yy_mb(i) - yy_mb(i-1))^2);
    end
end
pc_mb3 = pc_mb3 + sqrt((xx_mb(1) - xx_mb(end))^2 +...
    (yy_mb(1) - yy_mb(end))^2);

plot([xx_mb;xx_mb(1)],[yy_mb;yy_mb(1)],'Color','r')
fprintf('pc_mb3 = %f\n', pc_mb3);

% [x1_mb, x2_mb] = ndgrid(x_mb, x_mb);
% [y1_mb, y2_mb] = ndgrid(y_mb, y_mb);
% d_mb = ((x1_mb - x2_mb).^2 + (y1_mb - y2_mb).^2).^(1/2);
% d_mb(d_mb == 0) = inf;
% ind_mb = zeros(n_mb,n_mb);
% for i = 1 : n_mb
%     ind_mb(i,:) = d_mb(i,:) == min(d_mb(i,:));
% end

% edg_skl = bwmorph(edg_dil, 'skel', Inf);
% pc_edgskl = nnz(edg_skl);
% fprintf('pc_edgskl = %f\n', pc_edgskl);

figure(2)
img_eroded = imerode(img,SE);
edg_erd = img_eroded - img;
pc_edgerd = nnz(edg_erd);
fprintf('pc_edgerd = %f\n', pc_edgerd);
% imshow(edg_erd);
imov2 = imoverlay(img, edg_erd, 'b');
imshow(imov2);

% pc_final = get_perimeter2(~img);
% fprintf('pc_final = %f\n', pc_final);

%== GET_PERIMETER2 =============================================================%
%   An updated method to get the perimeter of the aggregate.
%   AUTHOR:  Hamed Nikookar, Timothy Sipkens, 2022-03-25
function p = get_perimeter2(img_binary)
 
mb = bwboundaries(img_binary);  % time-limiting step
 
x_mb = mb{1}(:,2);  % get x and y coordinates
y_mb = mb{1}(:,1);
 
fx = [0; x_mb(2:end) ~= x_mb(1:end-1)];  % flag change in x
fy = [0; y_mb(2:end) ~= y_mb(1:end-1)];  % flag change in y
ii_mb = cumsum(and(fx, fy)) + 1;  % count diagonal elements
 
% Combine first and last elements
if (x_mb(1) == x_mb(end)) || (y_mb(1) == y_mb(end))
    ii_mb(ii_mb == ii_mb(end)) = 1;
end
 
% Average x and y over each line segment.
n2 = accumarray(ii_mb, ones(length(ii_mb)));
x2 = accumarray(ii_mb, x_mb) ./ n2;
y2 = accumarray(ii_mb, y_mb) ./ n2;
p = sum(sqrt(diff(x2) .^ 2 + diff(y2) .^ 2));
p = p + ...
    sqrt((x2(end) - x2(1)) .^ 2 + ...
    (y2(end) - y2(1)) .^ 2);  % join outline
 
end
