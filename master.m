function [op3, tx, ty] = master(im, file)
%MASTER master funtion for flow based cartoonification
%   takes as imput a color image (RGB) im
%   returns output image after cartoonification (op3) 
%   along with the ETF tx, ty for debugging

[tx, ty] = DIP_ETF(im, 0);
% [tx2, ty2] = DIP_ETF(im, 2);
% [tx3, ty3] = DIP_ETF(im, 3);

op1 = lineExtract(tx, ty, rgb2gray(im));
% old code
% o2 = lineExtract(tx1, ty1, im(:,:,2));
% o3 = lineExtract(tx1, ty1, im(:,:,3));
% op3 = o1.*o2.*o3;s
% [tx, ty] = DIP_ETF(im, 0);

op2 = colorSmoothen(im, tx,  ty);

op3 = op1.*op2;

imwrite(op3, file);

end

