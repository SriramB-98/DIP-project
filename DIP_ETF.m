function [TdirX, TdirY] = DIP_ETF(im, ch)

% im = rgb2ycbcr(im);
% im = rgb2gray(im);
if ch > 0
    im = im(:, :, ch);
else
    im = rgb2gray(im);
end
img = im2double(im);
%rgb2gray(im));

[r, c, ~] = size(img);
A = imgaussfilt(img, 2);

mu = 5;
n = 3;

[Gx, Gy] = imgradientxy(A);
[Gmag, ~] = imgradient(Gx, Gy);
TdirX = Gy./Gmag;
% cos(Tdir);
TdirY = -Gx./Gmag;
% sin(Tdir);
Gmag = Gmag/max(max(Gmag));

for a = 1:n

    Tdir_newX = zeros(r, c);
    Tdir_newY = zeros(r, c);

    for i = 1:r
        for j = 1:c
            il = max(i - mu, 1);
            iu = min(i + mu, r);
            jl = max(j - mu, 1);
            ju = min(j + mu, c);
			
            counter = 0;
            
            for it = il:iu
                for jt = jl:ju
%                     if((it - i)*(it - i) + (jt - j)*(jt - j) < mu*mu)
                        coeff = (TdirX(it,jt)*TdirX(i,j) + TdirY(it,jt)*TdirY(i,j))*(Gmag(it,jt) - Gmag(i,j) + 1)/2;
                        Tdir_newX(i,j) = Tdir_newX(i,j) + coeff*TdirX(it,jt);
                        Tdir_newY(i,j) = Tdir_newY(i,j) + coeff*TdirY(it,jt);
                        counter = counter + 1;
%                     end
                end
            end
           
        end
    end

    Tmag = sqrt(Tdir_newX.*Tdir_newX + Tdir_newY.*Tdir_newY);
    %maxmag = max(max(Tmag));
    TdirX = Tdir_newX./Tmag;
    TdirY = Tdir_newY./Tmag;
end

