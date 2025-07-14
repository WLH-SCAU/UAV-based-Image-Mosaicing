%   Distribution code Version 1.0 -- 06/31/2025 by Linhui Wang Copyright 2025
%
%   The code is created based on the method described in the following paper 
%   [1] ""UAV-based Image Mosaicing system for Agricultural Applications Using Novel B-SIFT-ILS Algorithm"*. Linhui Wang, Yongda Lin, Zhenqi Zhou, Xuxiang Peng, Lizhi Chen, Quanli Tang, and Yonghong Tan, IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 
%        presented at 2025. 
%  
%   The code and the algorithm are for non-comercial use only.

% Before running, please ensure that VLFeat is installed and set up properly
run('D:\Matlab\vlfeat-0.9.21\toolbox\vl_setup'); %Please change to your toolbox directory

%----------------------------
clear;
clc;
close all;

img1 = imread('001_IMG.JPG');
img2 = imread('002_IMG.JPG');
% Read two images

% Convert to grayscale
gray1 = single(rgb2gray(img1));
gray2 = single(rgb2gray(img2));

% Extract SIFT features
[f1, d1] = vl_sift(gray1);
[f2, d2] = vl_sift(gray2);

% Feature matching
[matches, scores] = vl_ubcmatch(d1, d2);

% Extract matched points
matchedPoints1 = f1(1:2, matches(1,:))';
matchedPoints2 = f2(1:2, matches(2,:))';

% Estimate homography matrix using estimateGeometricTransform2D + RANSAC (i.e., BANSAC)
[tform, inlierIdx] = estimateGeometricTransform2D(matchedPoints1, matchedPoints2, ...
    'projective', 'MaxNumTrials', 2000, 'Confidence', 99.9, 'MaxDistance', 4);

% Extract inlier matches
inlierPoints1 = matchedPoints1(inlierIdx, :);
inlierPoints2 = matchedPoints2(inlierIdx, :);

% Display matching results (only inliers)
figure;
showMatchedFeatures(img1, img2, inlierPoints1, inlierPoints2, 'montage');
title('BANSAC Inlier Matches');

% Display the estimated homography matrix
disp('Estimated Homography Matrix:');
disp(tform.T);

% Image warping
outputView = imref2d(size(img2));
warpedImage1 = imwarp(img1, tform, 'OutputView', outputView);

% Create stitched image
blended = max(warpedImage1, img2); % Simple blending method

% Display the stitching result
figure;
imshow(blended);
title('Stitched Result');
