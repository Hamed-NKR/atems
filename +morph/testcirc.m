clc
clear
close all

n_img = [1024, 1024];
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

figure(2)
img_eroded = imerode(img,SE);
edg_erd = img_eroded - img;
pc_edgerd = nnz(edg_erd);
fprintf('pc_edgerd = %f\n', pc_edgerd);
% imshow(edg_erd);
imov2 = imoverlay(img, edg_erd, 'b');
imshow(imov2);

% edg_skl = bwmorph(edg_dil, 'skel', Inf);
% pc_edgskl = nnz(edg_skl);
% fprintf('pc_edgskl = %f\n', pc_edgskl);

mb = bwboundaries(~img);
pc_mb1 = size(mb{1},1);
fprintf('pc_mb1 = %f\n', pc_mb1);

n_mb = length(mb{1});
pc_mb2 = 0;
for i = 1 : n_mb
    if i ~= n_mb
        pc_mb2 = pc_mb2 + sqrt(sum((mb{1}(i+1,:) - mb{1}(i,:)).^2));
    else
        pc_mb2 = pc_mb2 + sqrt(sum((mb{1}(1,:) - mb{1}(i,:)).^2));
    end
end
fprintf('pc_mb2 = %f\n', pc_mb2);
    
    
    
    
