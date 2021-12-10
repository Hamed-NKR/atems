% BLUR  Blurs a transmission map image by the given radius w. 
%
% INPUTS:  
%   [A] - Image array to blur
%   [r] - Blur radius
%
% OUTPUTS:
%   [OUTPUT] - Blurred image array
% Darwin Zhu, 2021-09-26
%=========================================================================%

function [output] = blur(A,r)
[row, col] = size(A);
B=ones(size(A) + (2*r));
B(r+1:end-r,r+1:end-r)=A;
output = 0*A;
for i=r+1:row+r
  for j=r+1:col+r
    tmp=B(i-r:i+r,j-r:j+r);
    output(i-r,j-r)=mean(tmp(~isnan(tmp)));
  end
end