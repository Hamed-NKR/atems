clc
clear
close all

n_img = [1024, 1024];
% pixsz = 1.9;

img = ones(n_img(1), n_img(2));

c = [randperm(n_img(1), 1), randperm(n_img(2), 1)];
r1 = min(min(n_img(1) - c(1), c(1) - 1),...
    min(n_img(2) - c(2), c(2) - 1)) / 2;
r2 = min(min(n_img(1) - c(1), c(1) - 1),...
    min(n_img(2) - c(2), c(2) - 1)) / 5;

for i = 1 : n_img(1)
    for j = 1 : n_img(2)
        
        if (i-c(1))^2 / r1^2 + (j-c(2))^2 / r2^2 <= 1
            img(i,j) = 0;
        end
    end
end

a_exact = pi * r1 * r2;
fprintf('a_exact = %f\n', a_exact);
a1_img = numel(img) - nnz(img);
fprintf('a1_img = %f\n', a1_img);
a2_img = bwarea(~img);
fprintf('a_img = %f\n', a2_img);

e = sqrt(r1^2 - r2^2) / r1;
h = (r1 - r2)^2 / (r1 + r2)^2;
p1_exact = 2 * r1 * pi * (1 - (1/2)^2 * e^2 - ((1*3)/(2*4))^2 * e^4 / 3 -...
    ((1*3*5)/(2*4*6))^2 * e^6 / 5);
p2_exact = pi * (r1 + r2) * (1 + (1/4) * h + (1/64) * h^2 + (1/256) * h^3);
fprintf('p1_exact = %f\n', p1_exact);
fprintf('p2_exact = %f\n', p2_exact);

figure(1)
imshow(img);

SE = strel('disk', 1);
img_dilated = imdilate(img,SE);
edg_dil = img_dilated - img;
pc_edgdil = nnz(edg_dil);
fprintf('pc_edgdil = %f\n', pc_edgdil);

figure(2)
imshow(edg_dil);
% imov1 = imoverlay(edg_dil, img, 'r');
% imshow(imov1);

img_eroded = imerode(img,SE);
edg_erd = img_eroded - img;
pc_edgerd = nnz(edg_erd);
fprintf('pc_edgerd = %f\n', pc_edgerd);

figure(3)
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

