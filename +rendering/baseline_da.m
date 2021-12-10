% BASELINE_DA  Returns the projected area of image by counting white pixels.
%  
% INPUTS:
%   [IMG] - Image matrix
%   [PIXSIZE] - Size of pixels in nanometers
%
% OUTPUTS:
%   [DA] - Aggregate projected area
%  
%  AUTHOR: Darwin Zhu, 2021-09-10
%=========================================================================%


function [da] = baseline_da(img, pixsize)
img_size = size(img);
x = img_size(1);
y = img_size(2);
count = x*y;
for i = 1:x
for j = 1:y
if img(i,j) == 255
count = count - 1;
end
end
end
da = ((count/pi)^.5)*2*pixsize;