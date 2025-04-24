function varargout = compare_edge_detectors(img)
% COMPARE_EDGE_DETECTORS 比较不同边缘检测算法的效果
% 输入:
%   img - 输入图像 (灰度图)，如果没有提供，将提示用户选择一个图像文件

% 如果没有提供输入参数，提示用户选择图像
if nargin < 1
    [filename, pathname] = uigetfile({'*.jpg;*.png;*.bmp;*.tif', '图像文件 (*.jpg, *.png, *.bmp, *.tif)'}, '选择一个图像文件');
    if isequal(filename, 0) || isequal(pathname, 0)
        disp('用户取消了操作');
        return;
    end
    img = imread(fullfile(pathname, filename));
    disp(['已加载图像: ', fullfile(pathname, filename)]);
end

% 确保图像是灰度图
if size(img, 3) == 3
    img = rgb2gray(img);
end

% 使用更温和的预处理参数，减少过度模糊
preprocessed_img = imgaussfilt(img, 1.2);  % 降低sigma值从2.5到1.2，减少模糊

% 使用自适应直方图均衡来增强对比度，同时保留更多细节
preprocessed_img = adapthisteq(preprocessed_img, 'ClipLimit', 0.02, 'Distribution', 'rayleigh');

% 使用较小的窗口进行中值滤波，保留边缘细节
preprocessed_img = medfilt2(preprocessed_img, [2 2]);  % 减小窗口大小从3x3到2x2

% Canny算法 - 使用适度的阈值
[~, thresh] = edge(preprocessed_img, 'Canny');
thresh_high = thresh(2) * 1.2;  % 降低高阈值乘数，从2.5改为1.2
thresh_low = thresh(1) * 0.8;   % 降低低阈值乘数，从2.0改为0.8
edges_canny = edge(preprocessed_img, 'Canny', [thresh_low thresh_high]);

% 其他边缘检测算法
edges_roberts = edge(img, 'Roberts');
edges_sobel = edge(img, 'Sobel');
edges_prewitt = edge(img, 'Prewitt');
edges_log = edge(img, 'log');

% 改进 Canny 边缘检测以增强内部圈孔的检测能力

% 1. 对图像进行局部对比度增强，突出圆孔边缘
preprocessed_img = adapthisteq(preprocessed_img, 'ClipLimit', 0.02, 'Distribution', 'rayleigh');

% 2. 使用温和的高斯滤波，平衡噪声抑制和边缘保留
preprocessed_img = imgaussfilt(preprocessed_img, 1.5);

% 3. 使用更激进的双阈值 Canny，确保捕获更多弱边缘
edges_canny = edge(preprocessed_img, 'Canny', [0.04 0.15]); % 非常低的阈值以确保捕获圆孔边缘

% 3. 形态学处理: 先闭运算连接断开的边缘
se1 = strel('disk', 1);
edges_closed = imclose(edges_canny, se1);

% 4. 再使用开运算去除小噪点
se2 = strel('disk', 2);
edges_cleaned = imopen(edges_closed, se2);

% 5. 在预处理图像上寻找圈形结构，优化半径范围和灵敏度
[centers, radii] = imfindcircles(preprocessed_img, [10 25], 'Sensitivity', 0.9, 'EdgeThreshold', 0.1);

% 6. 区域特异性处理：左右区域分别处理，右侧更严格筛选
[height, width] = size(edges_cleaned);
left_region = edges_cleaned(:, 1:round(width/2));
right_region = edges_cleaned(:, round(width/2)+1:end);

% 左区域使用小面积阈值
left_region = bwareaopen(left_region, 20);

% 右区域使用更大面积阈值和形状分析
right_region = bwareaopen(right_region, 60);

% 对右区域进行连通区域分析，筛选形状
cc_right = bwconncomp(right_region);
stats_right = regionprops(cc_right, 'Area', 'Perimeter', 'Solidity');

% 设置筛选条件：移除面积小且周长/面积比大的区域（非圆形）
for i = 1:length(stats_right)
    if stats_right(i).Area < 70 && (stats_right(i).Perimeter / stats_right(i).Area > 0.7)
        % 移除不符合条件的区域
        right_region(cc_right.PixelIdxList{i}) = 0;
    end
end

% 合并处理后的左右区域
edges_canny = [left_region, right_region];

% 7. 如果找到圈形，将它们叠加到edges_canny中增强可见性
if ~isempty(centers)
    disp(['检测到 ', num2str(length(radii)), ' 个圈形结构']);
    circle_mask = false(size(preprocessed_img));
    for i = 1:length(radii)
        [x, y] = circlepoints(round(centers(i,1)), round(centers(i,2)), round(radii(i)));
        x = max(1, min(size(circle_mask, 2), x));
        y = max(1, min(size(circle_mask, 1), y));
        indices = sub2ind(size(circle_mask), y, x);
        circle_mask(indices) = true;
    end
    % 将圈形边缘叠加到结果中
    edges_canny = edges_canny | circle_mask;
end

% 打印检测到的像素数
disp(['原始Canny检测结果像素数：', num2str(sum(edges_canny(:)))]);

% 打印当前像素数
disp(['形态学处理后的Canny像素数：', num2str(sum(edges_canny(:)))]);

% 如果像素数未在理想范围内，进行调整
if sum(edges_canny(:)) < 800
    % 如果像素太少，降低阈值并尝试更强的闭合操作
    edges_canny = edge(preprocessed_img, 'Canny', [0.08 0.15]);
    se = strel('disk', 3);
    edges_canny = imclose(edges_canny, se);
    disp(['使用更强的处理后的Canny像素数：', num2str(sum(edges_canny(:)))]);
elseif sum(edges_canny(:)) > 3000
    % 如果像素太多，提高阈值
    edges_canny = edge(preprocessed_img, 'Canny', [0.15 0.3]);
    disp(['提高阈值后的Canny像素数：', num2str(sum(edges_canny(:)))]);
end

% 改进的Canny边缘处理，重点处理区域特异性噪声
[m, n] = size(edges_canny);

% 创建左右区域分隔点
right_start = round(n * 0.55); % 从图像55%处开始定义为“右侧”

% 1. 清除边缘像素
margin = 3; % 稍微增加边距
edges_canny_cleaned = edges_canny;
edges_canny_cleaned(1:margin, :) = false;
edges_canny_cleaned(m-margin+1:m, :) = false;
edges_canny_cleaned(:, 1:margin) = false;
edges_canny_cleaned(:, n-margin+1:n) = false;

% 2. 分区域处理
% 左侧区域 - 适度清理
left_region = edges_canny_cleaned(:, 1:right_start);
left_cleaned = bwareaopen(left_region, 20); % 降低阈值以保留小结构

% 右侧区域 - 适度清理，但要注意保留圈形结构
right_region = edges_canny_cleaned(:, right_start+1:end);
% 先检测有没有圈形结构
cc_right = bwconncomp(right_region);
stats_right = regionprops(cc_right, 'Area', 'Perimeter', 'Extent', 'Solidity');

% 找出可能的圈形结构（面积中等，密集度高）
circle_mask = false(size(right_region));
for i = 1:length(stats_right)
    circularity = 4 * pi * stats_right(i).Area / ((stats_right(i).Perimeter)^2);
    % 圈形的特征: 圈形度高，密度高
    if circularity > 0.6 && stats_right(i).Solidity > 0.7
        circle_mask(cc_right.PixelIdxList{i}) = true;
    end
end

% 对非圈形结构进行面积过滤
right_non_circle = right_region & ~circle_mask;
right_non_circle_filtered = bwareaopen(right_non_circle, 50);

% 合并圈形和其他过滤结构
right_cleaned = right_non_circle_filtered | circle_mask;

% 输出圈形检测信息
disp(['检测到的可能圈形结构像素数：', num2str(sum(circle_mask(:)))]);

% 合并区域
edges_canny_cleaned = zeros(size(edges_canny_cleaned), 'logical');
edges_canny_cleaned(:, 1:right_start) = left_cleaned;
edges_canny_cleaned(:, right_start+1:end) = right_cleaned;

% 4. 对左右区域进行不同的形态学处理
% 左侧主要进行边缘闭合和连接
se_left = strel('disk', 2);
left_region = edges_canny_cleaned(:, 1:right_start);
left_processed = imclose(left_region, se_left);

% 右侧使用更复杂的处理，主要是清除杂象同时保留圈孔
% 先进行形态学给小线段连接的机会
right_region = edges_canny_cleaned(:, right_start+1:end);
se_right = strel('disk', 1);
right_processed = imclose(right_region, se_right);

% 对右侧区域进行连通组件分析，区分处理不同形状
cc = bwconncomp(right_processed);
stats = regionprops(cc, 'Area', 'Perimeter', 'Circularity');

% 形状分析，区分处理线段和圈形
% 1. 计算圈形度 = 4*pi*Area/(Perimeter^2)
% 2. 圈形度接近1的是圈形，微小的是线段
for i = 1:length(stats)
    % 计算圈形度（如枟 regionprops 不支持 Circularity）
    if ~isfield(stats, 'Circularity')
        circularity = 4 * pi * stats(i).Area / (stats(i).Perimeter^2);
    else
        circularity = stats(i).Circularity;
    end
    
    % 只过滤细长线段，保留圈形结构（圈形度>0.6）
    if stats(i).Area < 70 && circularity < 0.6 && (stats(i).Perimeter / stats(i).Area) > 0.7
        right_processed(cc.PixelIdxList{i}) = false;
    end
end

% 合并左右区域
edges_canny_final = false(size(edges_canny));
edges_canny_final(:, 1:right_start) = left_processed;
edges_canny_final(:, right_start+1:end) = right_processed;

% 5. 最终面积过滤保证像素数量合适
% 将阈值降低到 20，以保留小圈孔
edges_canny_final = bwareaopen(edges_canny_final, 20);

% 打印处理过程中的像素数量变化
disp(['初始的Canny像素数：', num2str(sum(edges_canny(:)))]);
disp(['区域清理后的像素数：', num2str(sum(edges_canny_cleaned(:)))]);
disp(['最终Canny像素数：', num2str(sum(edges_canny_final(:)))]);

% 确保像素数量在合理范围内 (1000-1500)
if sum(edges_canny_final(:)) > 1600
    % 递增面积阈值直到像素数在合适范围
    area_threshold = 40;
    while sum(edges_canny_final(:)) > 1400 && area_threshold < 100
        area_threshold = area_threshold + 10;
        edges_canny_final = bwareaopen(edges_canny_final, area_threshold);
    end
    disp(['面积调整后的像素数（阈值=', num2str(area_threshold), '）：', ...
          num2str(sum(edges_canny_final(:)))]);
end

% 所有结果显示在单一窗口中 - 3x3网格
figure('Name', '边缘检测算法全面比较', 'Position', [50, 50, 1200, 800]);

% 打印所有算法的像素计数信息
disp(['Canny结果像素数：', num2str(sum(edges_canny_final(:)))]);
disp(['Roberts结果像素数：', num2str(sum(edges_roberts(:)))]);
disp(['Sobel结果像素数：', num2str(sum(edges_sobel(:)))]);
disp(['Prewitt结果像素数：', num2str(sum(edges_prewitt(:)))]);
disp(['LoG结果像素数：', num2str(sum(edges_log(:)))]);

% 第一行：原始图像、预处理图像、Canny结果
subplot(3, 3, 1);
imshow(img);
title('原始图像');

subplot(3, 3, 2);
imshow(preprocessed_img);
title('预处理后的图像');

subplot(3, 3, 3);
imshow(edges_canny_final);
title(['Canny: ', num2str(sum(edges_canny_final(:))), ' 像素']);

% 第二行：其他边缘检测算法结果
subplot(3, 3, 4);
imshow(edges_roberts);
title(['Roberts: ', num2str(sum(edges_roberts(:))), ' 像素']);

subplot(3, 3, 5);
imshow(edges_sobel);
title(['Sobel: ', num2str(sum(edges_sobel(:))), ' 像素']);

subplot(3, 3, 6);
imshow(edges_prewitt);
title(['Prewitt: ', num2str(sum(edges_prewitt(:))), ' 像素']);

% 第三行：LoG结果和叠加效果
subplot(3, 3, 7);
imshow(edges_log);
title(['LoG: ', num2str(sum(edges_log(:))), ' 像素']);

% Canny结果叠加在原图上
subplot(3, 3, 8);
imshow(img);
hold on;
[edge_y, edge_x] = find(edges_canny_final);
plot(edge_x, edge_y, 'r.', 'MarkerSize', 1);
hold off;
title('Canny结果叠加在原图上');

% LoG结果叠加在原图上
subplot(3, 3, 9);
imshow(img);
hold on;
[log_y, log_x] = find(edges_log);
plot(log_x, log_y, 'g.', 'MarkerSize', 1);
hold off;
title('LoG结果叠加在原图上');

% 返回不同算法的结果供进一步分析
if nargout > 0
    varargout{1} = edges_canny_final;  % 单一Canny结果
    varargout{2} = edges_roberts;
    varargout{3} = edges_sobel;
    varargout{4} = edges_prewitt;
    varargout{5} = edges_log;
end
end
