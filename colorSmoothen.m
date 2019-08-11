function [o_image] = colorSmoothen(img, TdirX, TdirY)
%COLORSMOOTHEN Implements Flow based Bilateral Filtering
%   This function takes in the ETF and caluculates the flow based bilateral
%   filtering by moving along and perpendicular to the ETF

% Discretises the angles
angle = atan2d(TdirY,TdirX);
angle = round(angle/45)*45;

[r, c, d] = size(img);
length = 4; %length of path to be traced

img = im2double(img);
o_image = img;
% sigma for spacial and intensity domain
sigma_s_t = 5;
sigma_i_t = 50;
sigma_s_p = 2;
sigma_i_p = 10;

for ch = 1:d
    img = o_image(:,:,ch);
for k = 1:5
    % loop along the tangent curve
    new_img = zeros(r,c);
    for i = 1:r
        for j = 1:c
            num = 0;
            den = 0;
            s=0;
            i_n = i;
            j_n = j;
            i_nn = 0;
            j_nn = 0;
                
            %loop over the lenght required
            for z = 1:length
                % updating the next pixel and the condition for no update
                % This is for 1 to length
                if angle(i_n,j_n)==0
                    i_nn = i_n+1;
                    j_nn = j_n;
                elseif abs(angle(i_n,j_n))==45
                        i_nn = i_n + 1;
                        j_nn = j_n - abs(angle(i_n,j_n))/angle(i_n,j_n);
                elseif abs(angle(i_n,j_n)) == 90
                        j_nn = j_n - abs(angle(i_n,j_n))/angle(i_n,j_n);
                        i_nn = i_n;
                elseif abs(angle(i_n,j_n))==135
                        i_nn = i_n - 1;
                        j_nn = j_n - abs(angle(i_n,j_n))/angle(i_n,j_n);
                elseif abs(angle(i_n,j_n)) == 180
                        i_nn = i_n - 1; 
                        j_nn = j_n;
                end
                if (i_nn>r || j_nn>c || i_nn<1 || j_nn<1)
                     break;
                else
                     s = s+1;
                     i_n = i_nn;
                     j_n = j_nn;
                end
                % Bilateral filtering function
                num = num + img(i_n,j_n,:)*exp(-(img(i_n,j_n)-img(i,j))/2/sigma_i_t^2)*exp(-((i-i_n)^2+(j-j_n)^2)/2/sigma_s_t^2);
                den = den + exp(-(img(i_n,j_n)-img(i,j))/2/sigma_i_t^2)*exp(-((i-i_n)^2+(j-j_n)^2)/2/sigma_s_t^2);
            end
            
            i_n=i;
            j_n=j;
            
            for z = 1:length
                % updating the next pixel and the condition for no update
                % This is for -length to -1
                if angle(i_n,j_n)==0
                    i_nn = i_n-1;
                    j_nn = j_n;
                elseif abs(angle(i_n,j_n))==45
                        i_nn = i_n - 1;
                        j_nn = j_n + abs(angle(i_n,j_n))/angle(i_n,j_n);
                elseif abs(angle(i_n,j_n)) == 90
                        j_nn = j_n + abs(angle(i_n,j_n))/angle(i_n,j_n);
                        i_nn = i_n;
                elseif abs(angle(i_n,j_n))==135
                        i_nn = i_n + 1;
                        j_nn = j_n + abs(angle(i_n,j_n))/angle(i_n,j_n);
                elseif abs(angle(i_n,j_n)) == 180
                        i_nn = i_n + 1; 
                        j_nn = j_n;
                end
                if (i_nn>r || j_nn>c || i_nn<1 || j_nn<1)
                     break;
                else
                     s=s+1;
                     i_n = i_nn;
                     j_n = j_nn;
                end
                num = num + img(i_n,j_n)*exp(-(img(i_n,j_n)-img(i,j))^2/2/sigma_i_t^2)*exp(-((i-i_n)^2+(j-j_n)^2)/2/sigma_s_t^2);
                den = den + exp(-(img(i_n,j_n)-img(i,j))^2/2/sigma_i_t^2)*exp(-((i-i_n)^2+(j-j_n)^2)/2/sigma_s_t^2);
            end
            
            % summation of the bilteral filter
            num=num+img(i,j);
            den=den+1;
            new_img(i,j) = num/den;
            
        end
    end

    %perpendicular angles  [-90,-45,0,45,90] will become [0,45,90,135,180]
    perp_angle = 90 + angle;
    perp_angle(perp_angle>180) = perp_angle(perp_angle>180)-180;
    new_img_p = zeros(r, c);

    %same loop as above but the intensities considered is wrt new_img
    for i = 1:r
        for j = 1:c
            num = 0;
            den = 0;
            s=0;
            i_n = i;
            j_n = j;
            for z = 1:length
                if perp_angle(i_n,j_n)==0
                    i_nn = i_n+1;
                    j_nn = j_n;
                elseif abs(perp_angle(i_n,j_n))==45
                        i_nn = i_n + 1;
                        j_nn = j_n - abs(perp_angle(i_n,j_n))/perp_angle(i_n,j_n);
                elseif abs(angle(i_n,j_n)) == 90
                        j_nn = j_n - abs(perp_angle(i_n,j_n))/perp_angle(i_n,j_n);
                        i_nn = i_n;
                elseif abs(angle(i_n,j_n))==135
                        i_nn = i_n - 1;
                        j_nn = j_n - abs(perp_angle(i_n,j_n))/perp_angle(i_n,j_n);
                elseif abs(angle(i_n,j_n)) == 180
                        i_nn = i_n - 1;
                        j_nn = j_n;
                end
                if (i_nn>r || j_nn>c || i_nn<1 || j_nn<1)
                     break;
                else
                     s=s+1;
                     i_n = i_nn;
                     j_n = j_nn;
                end
                num = num + new_img(i_n,j_n)*exp(-(new_img(i,j)-new_img(i_n,j_n))^2/2/sigma_i_p^2)*exp(-((i-i_n)^2+(j-j_n)^2)/2/sigma_s_p^2);
                den = den + exp(-(new_img(i,j)-new_img(i_n,j_n))^2/2/sigma_i_p^2)*exp(-((i-i_n)^2+(j-j_n)^2)/2/sigma_s_p^2);
            end
            
            i_n=i;
            j_n=j;
            
            for z = 1:length
                if perp_angle(i_n,j_n)==0
                    i_nn = i_n-1;
                    j_nn = j_n;
                elseif abs(perp_angle(i_n,j_n))==45
                        i_nn = i_n - 1;
                        j_nn = j_n + abs(perp_angle(i_n,j_n))/perp_angle(i_n,j_n);
                elseif abs(angle(i_n,j_n)) == 90
                        j_nn = j_n + abs(perp_angle(i_n,j_n))/perp_angle(i_n,j_n);
                        i_nn =i_n;
                elseif abs(angle(i_n,j_n))==135
                        i_nn = i_n + 1;
                        j_nn = j_n + abs(perp_angle(i_n,j_n))/perp_angle(i_n,j_n);
                elseif abs(angle(i_n,j_n)) == 180
                        i_nn = i_n + 1;
                        j_nn = j_n;
                end
                if (i_nn>r || j_nn>c || i_nn<1 || j_nn<1)
                     break;
                else
                     s=s+1;
                     i_n = i_nn;
                     j_n = j_nn;
                end
                num = num + new_img(i_n,j_n)*exp(-(new_img(i,j)-new_img(i_n,j_n))^2/2/sigma_i_p^2)*exp(-((i-i_n)^2+(j-j_n)^2)/2/sigma_s_p^2);
                den = den + exp(-(new_img(i,j)-new_img(i_n,j_n))^2/2/sigma_i_p^2)*exp(-((i-i_n)^2+(j-j_n)^2)/2/sigma_s_p^2);
            end
            num=num+img(i,j);
            den=den+1;
            new_img_p(i,j) = num/den;
        end
    end
    img = new_img_p;
end
o_image(:, :, ch) = img;
end
