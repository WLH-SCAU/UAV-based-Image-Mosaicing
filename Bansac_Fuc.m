%   Distribution code Version 1.0 -- 06/31/2025 by Linhui Wang Copyright 2025
%
%   The code is created based on the method described in the following paper 
%   [1] ""UAV-based Image Mosaicing system for Agricultural Applications Using Novel B-SIFT-ILS Algorithm"*. Linhui Wang, Yongda Lin, Zhenqi Zhou, Xuxiang Peng, Lizhi Chen, Quanli Tang, and Yonghong Tan, IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 
%        presented at 2025. 
%  
%   The code and the algorithm are for non-comercial use only.

function [bestH, bestInliers] = bansacH(X1, X2, N, threshold)
    % X1, X2: 2xN point coordinates
    % N: number of iterations
    % threshold: inlier threshold

    numPoints = size(X1, 2);
    bestInliers = [];
    bestH = eye(3);

    for i = 1:N
        % Randomly select 4 matching points
        idx = randperm(numPoints, 4);
        pts1 = X1(:, idx);
        pts2 = X2(:, idx);

        % Estimate homography matrix
        H = computeH(pts1, pts2);

        % Transform all points
        projX2 = H * [X2; ones(1, numPoints)];
        projX2 = projX2 ./ projX2(3, :);

        % Compute error
        errors = sqrt(sum((X1 - projX2(1:2, :)).^2, 1));

        % Select inliers
        inliers = find(errors < threshold);

        if length(inliers) > length(bestInliers)
            bestInliers = inliers;
            bestH = H;
        end
    end
end

function H = computeH(pts1, pts2)
    % Estimate homography matrix (DLT)
    n = size(pts1, 2);
    A = zeros(2*n, 9);

    for i = 1:n
        x = pts2(1,i); y = pts2(2,i);
        xp = pts1(1,i); yp = pts1(2,i);
        A(2*i-1,:) = [-x, -y, -1, 0, 0, 0, x*xp, y*xp, xp];
        A(2*i,:) = [0, 0, 0, -x, -y, -1, x*yp, y*yp, yp];
    end

    [~, ~, V] = svd(A);
    H = reshape(V(:,9), 3, 3)';
end
