function [image_out] = lineExtract(Tx, Ty, image)
%LINEEXTRACT Extracts lines from the ETF Tx,Ty
% expects a grayscale image

% initialize size and convert ot double precis
image = im2double(image);
[r, c] = size(image);

% set parameters
niter = 3;
sigm = 3;
sigc = 1;
sigs = 1.6*sigc;
ro = 0.987;
tau = 0.7;
Smax = 5;
% ceil(sigm*4);
Tmax = 4;
% ceil(sigc*3);

% Gx = -Ty;
% Gy = Tx;

% extract gradient from image
[Gx, Gy] = imgradientxy(image);

% old code --- ignore
% [~, Gdir] = imgradient(Gx, Gy);
% mu = Tmax;
% img = image;
% image = zeros(r+2*mu, c+2*mu);
% image(mu+1:mu+r, mu+1:mu+c) = img;


% iterative construction of filter H
for iter = 1:niter
	
	% initialize H_g
    H_g = zeros(r, c);

	% direction of gradient
    angle = atan2d(Gy,Gx);
    
	% angle is discretized to multiples of 45
    angle = round(angle/45)*45;

	% iterate over all points in image
    for i = 1:r
        for j = 1:c            
            trange = -Tmax:Tmax;
            gkern = normpdf(trange, 0, sigc) - ro*normpdf(trange, 0, sigs);
            count = 0;
            
            % taking possible values of angle case by case
            if(angle(i,j) == 0 || abs(angle(i,j) == 180))
                for t = -Tmax:Tmax
                    if(i+t >= 1 && i+t <= r) 
                        H_g(i,j) = H_g(i,j) + image(i+t, j)*gkern(t+Tmax+1);
                        count = count + gkern(t+Tmax+1);
                    end
                end
            elseif(abs(angle(i,j)) == 45)
                for t = -Tmax:Tmax
                    if(i+t >= 1 && i+t <= r && j+t >= 1 && j+t <= c) 
                        H_g(i,j) = H_g(i,j) + image(i+t, j+t)*gkern(t+Tmax+1);
                        count = count + gkern(t+Tmax+1);
                    end
                end
            elseif(abs(angle(i,j)) == 90)
                for t = -Tmax:Tmax
                    if(j+t >= 1 && j+t <= c) 
                        H_g(i,j) = H_g(i,j) + image(i, j+t)*gkern(t+Tmax+1);
                        count = count + gkern(t+Tmax+1);
                    end
                end
            elseif (abs(angle(i,j)) == 135)
                for t = -Tmax:Tmax
                    if(i-t <= r && i-t >= 1 && j+t >= 1 && j+t <= c) 
                        H_g(i,j) = H_g(i,j) + image(i-t, j+t)*gkern(t+Tmax+1);
                        count = count + gkern(t+Tmax+1);
                    end
                end
            end
            
            % setiing value of H_g (inner integral in equation(6))
            H_g(i,j) = H_g(i,j)/count;
        end
    end
    
    % initialize H_e
    H_e = zeros(r, c);
    
    % direction of tangent (flow)
    angle = atan2d(Ty, Tx);
    angle = round(angle/45)*45;

	% iterate over all points in image
    for i = 1:r
        for j = 1:c
            
            % value at (i,j)
            acc = H_g(i,j);
            count = 1;
            
            % dir to proceed in forward and backward directions along the flow starting at (i,j)
            for dir = [-1 1]
                
                uptr = i;
                vptr = j;
            
                for s = 1:Smax
	            	% cases for the angle
                    if(angle(uptr,vptr) == 0)
                        uptr = uptr + 1*dir;
                    elseif (angle(uptr,vptr) == 180 || angle(uptr,vptr) == -180)
                        uptr = uptr - 1*dir;
                    elseif(angle(uptr, vptr) == 45)
                        uptr = uptr + 1*dir;
                        vptr = vptr + 1*dir;
                    elseif(angle(i,j) == 90)
                        vptr = vptr + 1*dir;
                    elseif (angle(i,j) == 135)
                        uptr = uptr - 1*dir;
                        vptr = vptr + 1*dir;
                    elseif(angle(uptr, vptr) == -45)
                        vptr = vptr - 1*dir;
                    elseif(abs(angle(i,j)) == -90)
                        vptr = vptr - 1*dir;
                    elseif (abs(angle(i,j)) == -135)
                        uptr = uptr - 1*dir;
                        vptr = vptr - 1*dir;
                    end
					
					% updating the uptr and vptr to next value along flow
                    uptr = min(r, max(1, uptr));
                    vptr = min(c, max(1, vptr));
					% update acc (numerator) and count (denominator)
                    acc = acc + H_g(uptr, vptr)*normpdf(s, 0, sigm);
                    count = count + normpdf(s, 0, sigm);
                end
            end
            
            % setting H_e (outer integral in equation (6))
            H_e(i,j) = acc/count;
        end
    end
    
    extlines = zeros(r, c);
    
    % generate bitmap
    for i = 1:r
        for j = 1:c
            if (H_e(i,j) < 0 && 1 + tanh(H_e(i,j)) < tau)
                extlines(i,j) = 0;
            else
                extlines(i,j) = 1;
            end
        end
    end
    
    % filter image with bitmap
    image = image.*extlines;

    image_out = extlines;
end

end

