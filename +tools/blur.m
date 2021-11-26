function [output] = blur(A,w)
[row, col] = size(A);
B=ones(size(A) + (2*w));
B(w+1:end-w,w+1:end-w)=A;
output = 0*A;
for i=w+1:row+w
  for j=w+1:col+w
    tmp=B(i-w:i+w,j-w:j+w);
    output(i-w,j-w)=mean(tmp(~isnan(tmp)));
  end
end