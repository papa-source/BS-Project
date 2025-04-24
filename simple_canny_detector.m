function [edges, pixel_count] = simple_canny_detector(img, sigma, low_thresh, high_thresh)
% SIMPLE_CANNY_DETECTOR 简单的Canny边缘检测器，专注于圆孔检测
%
% 输入参数:
%   img - 输入图像 (灰度图)，如果没有提供，将提示用户选择一个图像文件
%   sigma - 高斯滤波的sigma值，控制模糊程度，默认为1.5
%   low_thresh - Canny低阈值，默认为0.04
%   high_thresh - Canny高阈值，默认为0.15
%
% 输出:
%   edges - 检测到的边缘二值图像
%   pixel_count - 边缘像素数量
%
% 用法示例:
%   [edges, count] = simple_canny_detector();    % 提示选择图像
%   [edges, count] = simple_canny_detector(img); % 使用提供的图像
%   [edges, count] = simple_canny_detector(img, 1.8, 0.03, 0.12); % 自定义参数

% 如果没有提供输入图像，提示用户选择一个
if nargin < 1
    [filename, pathname] = uigetfile({'*.jpg;*.png;*.bmp;*.tif', '图像文件 (*.jpg, *.png, *.bmp, *.tif)'});
    if isequal(filename, 0) || isequal(pathname, 0)
        disp('用户取消了操作');
        edges = [];
        pixel_count = 0;
        return;
    end
    img = imread(fullfile(pathname, filename));
    disp(['已加载图像: ', fullfile(pathname, filename)]);
end

% 设置默认参数
if nargin < 4
    high_thresh = 0.15;
end
if nargin < 3
    low_thresh = 0.04;
end
if nargin < 2
    sigma = 1.5;
end

% 确保图像是灰度图
if size(img, 3) == 3
    img = rgb2gray(img);
end

% 1. 对图像进行预处理
% 局部对比度增强（自适应直方图均衡化）
enhanced_img = adapthisteq(img, 'ClipLimit', 0.02, 'Distribution', 'rayleigh');

% 高斯滤波平滑图像，减少噪声
blurred_img = imgaussfilt(enhanced_img, sigma);

% 2. 应用Canny边缘检测
edges = edge(blurred_img, 'Canny', [low_thresh high_thresh]);

% 3. 应用形态学操作以改进边缘
% 闭运算连接断开的边缘
se_close = strel('disk', 1);
edges = imclose(edges, se_close);

% 计算边缘像素数量
pixel_count = sum(edges(:));

% 显示结果
figure('Name', '简单Canny边缘检测器结果');

% 原始图像
subplot(2, 2, 1);
imshow(img);
title('原始图像');

% 预处理后的图像
subplot(2, 2, 2);
imshow(blurred_img);
title(['预处理图像 (sigma = ', num2str(sigma), ')']);

% Canny边缘检测结果
subplot(2, 2, 3);
imshow(edges);
title(['Canny: ', num2str(pixel_count), ' 像素, 阈值 = [', num2str(low_thresh), ' ', num2str(high_thresh), ']']);

% 叠加在原图上的结果
subplot(2, 2, 4);
imshow(img);
hold on;
[rows, cols] = find(edges);
plot(cols, rows, 'r.', 'MarkerSize', 1);
title('边缘叠加在原图上');

% 返回结果
end
