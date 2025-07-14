%   Distribution code Version 1.0 -- 06/31/2025 by Linhui Wang Copyright 2025
%
%   The code is created based on the method described in the following paper 
%   [1] ""UAV-based Image Mosaicing system for Agricultural Applications Using Novel B-SIFT-ILS Algorithm"*. Linhui Wang, Yongda Lin, Zhenqi Zhou, Xuxiang Peng, Lizhi Chen, Quanli Tang, and Yonghong Tan, IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 
%        presented at 2025. 
%  
%   The code and the algorithm are for non-comercial use only.


function U =ILS_LNorm(F, lambda, p, eps, iter)

F = single(F); % 'single' precision is very important to reduce the computational cost

gamma = 0.5 * p - 1;
c =  p * eps^gamma;

[N, M, D] = size(F);
sizeI2D = [N, M];

otfFx = psf2otf_Dx(sizeI2D); % equal to otfFx = psf2otf(fx, sizeI2D) where fx = [1, -1];
otfFy = psf2otf_Dy(sizeI2D); % equal to otfFy = psf2otf(fy, sizeI2D) where fy = [1; -1];

Denormin = abs(otfFx).^2 + abs(otfFy ).^2;
Denormin = repmat(Denormin, [1, 1, D]);
Denormin = 1 + 0.5 * c * lambda * Denormin;

U = F; % smoothed image

Normin1 = fft2(U);

for k = 1: iter
    
    % Intermediate variables \mu update, in x-axis and y-axis direction
    u_h = [diff(U,1,2), U(:,1,:) - U(:,end,:)];
    u_v = [diff(U,1,1); U(1,:,:) - U(end,:,:)];
        
    mu_h = c .* u_h - p .* u_h .* (u_h .* u_h + eps) .^ gamma;
    mu_v = c .* u_v - p .* u_v .* (u_v .* u_v + eps) .^ gamma;
    
    % Update the smoothed image U
    Normin2_h = [mu_h(:,end,:) - mu_h(:, 1,:), - diff(mu_h,1,2)];
    Normin2_v = [mu_v(end,:,:) - mu_v(1, :,:); - diff(mu_v,1,1)];
    
    FU = (Normin1 + 0.5 * lambda * (fft2(Normin2_h + Normin2_v))) ./ Denormin;
    U = real(ifft2(FU));
    
    Normin1 = FU;  % This helps to further enlarge the smoothing strength
    
end
