function hole_detection_system()
% 板材冲孔质量检测系统
% 包含图像预处理、边缘提取、圆孔识别与拟合、轮廓线识别与拟合四个主要模块

% 读取图像
[filename, pathname] = uigetfile({'*.jpg;*.png;*.bmp', '图像文件 (*.jpg, *.png, *.bmp)'});
if filename == 0
    return;
end
img = imread(fullfile(pathname, filename));

%% 1. 图像预处理模块
% 1.1 图像灰度化盘
if size(img, 3) == 3
    gray_img = rgb2gray(img);
else
    gray_img = img;
end

% 1.2 自适应中值滤波去噪
filtered_img = adaptive_median_filter(gray_img);

% 1.3 对比度增强
filtered_img = imadjust(filtered_img, stretchlim(filtered_img, [0.05 0.95]));

%% 2. 图像分割与ROI提取模块
% 2.1 使用多级阈值分割
T1 = graythresh(filtered_img);  % Otsu方法计算全局阈值
binary_img = imbinarize(filtered_img, T1);

% 确保目标为白色（值为1），背景为黑色（值为0）
if mean(binary_img(:)) < 0.5
    binary_img = ~binary_img;
end

% 2.2 形态学处理改善边缘
se1 = strel('disk', 2);  % 减小结构元素大小
% 先进行开运算去除噪点
binary_img = imopen(binary_img, se1);
% 再进行闭运算连接边缘
binary_img = imclose(binary_img, se1);

% 2.3 保留内部圆孔
% 标记连通区域
[L, num] = bwlabel(~binary_img, 8);
stats = regionprops(L, 'Area', 'Circularity');

% 找出圆形区域（圆孔）
is_hole = false(num, 1);
for i = 1:num
    if stats(i).Area > 30 && stats(i).Area < 1000 && stats(i).Circularity > 0.85  % 调整面积范围和圆度阈值
        is_hole(i) = true;
    end
end

% 重建二值图像，保留圆孔
binary_img_with_holes = binary_img;
for i = 1:num
    if is_hole(i)
        binary_img_with_holes(L == i) = 0;  % 将圆孔区域设为黑色
    end
end

% 2.4 ROI提取
% 移除小面积区域
binary_img = binary_img_with_holes;
binary_img = bwareaopen(binary_img, 50);  % 降低面积阈值以保留更多细节
roi = extract_roi(binary_img);

%% 3. 边缘提取模块
% 3.1 使用二值图像直接提取边缘
edges = bwperim(~binary_img_with_holes);  % 使用带孔的二值图像直接提取边缘

% 3.2 优化边缘
se = strel('disk', 1);
edges = imdilate(edges, se);  % 轻微膨胀以增强边缘
edges = bwmorph(edges, 'thin', Inf);  % 细化边缘
edges = bwareaopen(edges, 15);  % 减小面积阈值，保留更多细节

% 显示中间结果用于调试
figure('Name', '图像处理中间结果', 'NumberTitle', 'off');
subplot(3,2,1); imshow(gray_img); title('灰度图像');
subplot(3,2,2); imshow(filtered_img); title('滤波和对比度增强');
subplot(3,2,3); imshow(binary_img); title('二值化结果');
subplot(3,2,4); imshow(roi); title('ROI提取结果');
subplot(3,2,5); 
imshow(edges);
title('边缘检测结果（包含圆孔）');

% 显示直方图
subplot(3,2,6); 
imhist(filtered_img);
title('图像直方图');
hold on;
y_limits = ylim;
plot([T1*255 T1*255], [0 y_limits(2)], 'r-', 'LineWidth', 2);
legend('直方图', '阈值');
hold off;

%% 4. 圆孔识别与拟合模块
% 使用边缘检测结果进行圆检测
[centers, radii] = imfindcircles(edges, [10 30], ...  % 缩小半径范围到10-30像素
    'ObjectPolarity', 'bright', ...
    'Sensitivity', 0.85, ...  % 降低灵敏度
    'EdgeThreshold', 0.2, ...  % 提高边缘阈值
    'Method', 'PhaseCode');

% 如果没有检测到圆，尝试使用二值掩码
if isempty(centers)
    holes_mask = ~binary_img_with_holes & binary_img;
    [centers, radii] = imfindcircles(holes_mask, [10 30], ...  % 同样缩小半径范围
        'Sensitivity', 0.85, ...
        'EdgeThreshold', 0.2);
end

%% 5. 轮廓线识别与拟合模块
% 使用工件的二值图像边缘进行轮廓检测
workpiece_edges = bwperim(binary_img_with_holes);

% 获取ROI区域
stats = regionprops(binary_img_with_holes, 'BoundingBox', 'Area');
if ~isempty(stats)
    % 选择最大面积的区域作为工件
    [~, max_idx] = max([stats.Area]);
    bbox = stats(max_idx).BoundingBox;
    roi_x = round(bbox(1));
    roi_y = round(bbox(2));
    roi_w = round(bbox(3));
    roi_h = round(bbox(4));
    
    % 创建ROI掩码，并添加边距
    margin = 10;  % 添加10像素的边距
    roi_mask = false(size(workpiece_edges));
    roi_y_min = max(1, roi_y - margin);
    roi_y_max = min(size(workpiece_edges,1), roi_y + roi_h + margin);
    roi_x_min = max(1, roi_x - margin);
    roi_x_max = min(size(workpiece_edges,2), roi_x + roi_w + margin);
    roi_mask(roi_y_min:roi_y_max, roi_x_min:roi_x_max) = true;
    
    % 只保留ROI区域内的边缘
    workpiece_edges = workpiece_edges & roi_mask;
    
    % 移除图像边框
    border_margin = 20;  % 设置边框区域的宽度
    workpiece_edges(1:border_margin,:) = 0;  % 清除上边框
    workpiece_edges(end-border_margin:end,:) = 0;  % 清除下边框
    workpiece_edges(:,1:border_margin) = 0;  % 清除左边框
    workpiece_edges(:,end-border_margin:end) = 0;  % 清除右边框
end

% 外轮廓检测
[H_outer, theta_outer, rho_outer] = hough(workpiece_edges);
P_outer = houghpeaks(H_outer, 12, 'threshold', ceil(0.25*max(H_outer(:))));  % 增加峰值点数量，降低阈值
lines = houghlines(workpiece_edges, theta_outer, rho_outer, P_outer, ...
    'FillGap', 25, 'MinLength', 25);  % 增加填充间隙，减小最小长度

% 合并相近的线段
if ~isempty(lines)
    merged_lines = [];
    used = false(1, length(lines));
    
    for i = 1:length(lines)
        if used(i)
            continue;
        end
        
        current_line = lines(i);
        used(i) = true;
        
        % 计算当前线段的角度
        dx = current_line.point2(1) - current_line.point1(1);
        dy = current_line.point2(2) - current_line.point1(2);
        angle1 = atan2(dy, dx);
        
        % 寻找可以合并的线段
        for j = i+1:length(lines)
            if used(j)
                continue;
            end
            
            % 计算待比较线段的角度
            dx = lines(j).point2(1) - lines(j).point1(1);
            dy = lines(j).point2(2) - lines(j).point1(2);
            angle2 = atan2(dy, dx);
            
            % 放宽合并条件
            angle_diff = abs(mod(angle1 - angle2 + pi, pi));
            if angle_diff < 0.3 && ...  % 增加角度容差
               (norm(current_line.point2 - lines(j).point1) < 40 || ...  % 增加距离容差
                norm(current_line.point1 - lines(j).point2) < 40)
                
                % 合并两条线段
                points = [current_line.point1; current_line.point2; 
                         lines(j).point1; lines(j).point2];
                [~, idx] = max(pdist2(mean(points), points));
                opposite_idx = mod(idx + 2, 4);
                if opposite_idx == 0
                    opposite_idx = 4;
                end
                
                current_line.point1 = points(idx, :);
                current_line.point2 = points(opposite_idx, :);
                used(j) = true;
            end
        end
        
        merged_lines = [merged_lines, current_line];
    end
    
    lines = merged_lines;
end

% 添加轮廓闭合处理
if ~isempty(lines)
    % 找到所有线段端点
    endpoints = [];
    for i = 1:length(lines)
        endpoints = [endpoints; lines(i).point1; lines(i).point2];
    end
    
    % 寻找未闭合的端点（与其他端点距离较远的点）
    n_endpoints = size(endpoints, 1);
    unclosed = false(n_endpoints, 1);
    for i = 1:n_endpoints
        min_dist = inf;
        for j = 1:n_endpoints
            if i ~= j
                dist = norm(endpoints(i,:) - endpoints(j,:));
                min_dist = min(min_dist, dist);
            end
        end
        if min_dist > 30  % 设置闭合阈值
            unclosed(i) = true;
        end
    end
    
    % 连接未闭合的端点
    unclosed_points = endpoints(unclosed,:);
    if size(unclosed_points, 1) >= 2
        for i = 1:2:size(unclosed_points, 1)
            if i+1 <= size(unclosed_points, 1)
                % 创建与原始线段结构相同的新线段
                new_line = struct('point1', unclosed_points(i,:), ...
                                'point2', unclosed_points(i+1,:), ...
                                'theta', 0, ...
                                'rho', 0);
                lines = [lines, new_line];
            end
        end
    end
end

% 获取角点
corners = [];
if length(lines) >= 2
    for i = 1:length(lines)-1
        for j = i+1:length(lines)
            corner = line_intersection(lines(i), lines(j));
            if ~isempty(corner)
                % 检查角点是否在ROI区域内
                if corner(1) >= roi_x && corner(1) <= roi_x+roi_w && ...
                   corner(2) >= roi_y && corner(2) <= roi_y+roi_h
                    corners = [corners; corner];
                end
            end
        end
    end
    
    % 合并相近的角点
    if ~isempty(corners)
        corners = merge_close_corners(corners);
    end
end

%% 6. 显示结果
display_results(img, edges, centers, radii, lines, corners);
end

function filtered = adaptive_median_filter(img)
% 自适应中值滤波实现
[m, n] = size(img);
filtered = zeros(m, n);
min_window = 3;
max_window = 7;

for i = 1:m
    for j = 1:n
        window_size = min_window;
        while window_size <= max_window
            window = get_window(img, i, j, window_size);
            z_min = min(window(:));
            z_max = max(window(:));
            z_med = median(window(:));
            z_xy = img(i,j);
            
            if z_med > z_min && z_med < z_max
                if z_xy > z_min && z_xy < z_max
                    filtered(i,j) = z_xy;
                else
                    filtered(i,j) = z_med;
                end
                break;
            else
                window_size = window_size + 2;
                if window_size > max_window
                    filtered(i,j) = z_med;
                end
            end
        end
    end
end
filtered = uint8(filtered);
end

function window = get_window(img, i, j, window_size)
% 获取图像窗口
half = floor(window_size/2);
[m, n] = size(img);
row_min = max(1, i-half);
row_max = min(m, i+half);
col_min = max(1, j-half);
col_max = min(n, j+half);
window = img(row_min:row_max, col_min:col_max);
end

function threshold = iterative_threshold(img)
% 迭代法求最佳阈值
threshold = mean(img(:));
delta = 0.1;
while true
    g1 = img(img >= threshold);
    g2 = img(img < threshold);
    m1 = mean(g1);
    m2 = mean(g2);
    new_threshold = (m1 + m2) / 2;
    if abs(new_threshold - threshold) < delta
        break;
    end
    threshold = new_threshold;
end
end

function roi = extract_roi(binary_img)
% ROI提取
stats = regionprops(binary_img, 'BoundingBox');
if ~isempty(stats)
    bbox = stats(1).BoundingBox;
    x = round(bbox(1));
    y = round(bbox(2));
    w = round(bbox(3));
    h = round(bbox(4));
    roi = binary_img(y:y+h-1, x:x+w-1);
else
    roi = binary_img;
end
end

function [centers, radii] = detect_and_fit_circles(edges)
% 圆孔识别与拟合
[centers, radii] = imfindcircles(edges, [20 100], 'Sensitivity', 0.85);
if isempty(centers)
    centers = [];
    radii = [];
end
end

function [lines, corners] = detect_and_fit_contours(edges)
% 轮廓线识别与拟合
[H, theta, rho] = hough(edges);
P = houghpeaks(H, 10, 'threshold', ceil(0.3*max(H(:))));
lines = houghlines(edges, theta, rho, P, 'FillGap', 5, 'MinLength', 7);

% 获取角点
corners = [];
if length(lines) >= 2
    for i = 1:length(lines)-1
        for j = i+1:length(lines)
            line1 = lines(i);
            line2 = lines(j);
            corner = line_intersection(line1, line2);
            if ~isempty(corner)
                corners = [corners; corner];
            end
        end
    end
end
end

function corner = line_intersection(line1, line2)
% 计算两直线交点
x1 = line1.point1(1);
y1 = line1.point1(2);
x2 = line1.point2(1);
y2 = line1.point2(2);
x3 = line2.point1(1);
y3 = line2.point1(2);
x4 = line2.point2(1);
y4 = line2.point2(2);

denominator = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4);
if abs(denominator) < 1e-10
    corner = [];
    return;
end

x = ((x1*y2-y1*x2)*(x3-x4) - (x1-x2)*(x3*y4-y3*x4)) / denominator;
y = ((x1*y2-y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4)) / denominator;
corner = [x, y];
end

function display_results(original_img, edges, centers, radii, lines, corners)
% 显示检测结果
figure('Name', '板材冲孔质量检测结果', 'NumberTitle', 'off');

% 原图
subplot(2,2,1);
imshow(original_img);
title('原始图像');

% 边缘检测结果
subplot(2,2,2);
imshow(edges);
title('边缘检测结果');

% 圆孔检测结果
subplot(2,2,3);
imshow(original_img);
hold on;
if ~isempty(centers)
    viscircles(centers, radii, 'EdgeColor', 'b', 'LineWidth', 1.5);
    % 标注圆心和半径
    for i = 1:size(centers, 1)
        plot(centers(i,1), centers(i,2), 'b+', 'MarkerSize', 10, 'LineWidth', 2);
        text(centers(i,1)+10, centers(i,2)+10, ...
             sprintf('R=%.1f', radii(i)), ...
             'Color', 'blue', 'FontSize', 8);
    end
end
title('圆孔检测结果');
hold off;

% 轮廓线检测结果
subplot(2,2,4);
imshow(original_img);
hold on;
% 绘制轮廓线
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1), xy(:,2), 'LineWidth', 2, 'Color', 'green');
end
% 绘制角点
if ~isempty(corners)
    plot(corners(:,1), corners(:,2), 'r*', 'MarkerSize', 10);
end
title('轮廓线检测结果');
hold off;
end

function merged_corners = merge_close_corners(corners)
% 合并相近的角点
if isempty(corners)
    merged_corners = corners;
    return;
end

% 设置距离阈值
dist_threshold = 10;

% 初始化标记数组
n = size(corners, 1);
merged = false(n, 1);
merged_corners = [];

% 遍历所有角点
for i = 1:n
    if merged(i)
        continue;
    end
    
    % 计算当前角点与其他角点的距离
    distances = sqrt(sum((corners - corners(i,:)).^2, 2));
    close_corners_idx = distances < dist_threshold;
    
    % 计算相近角点的平均位置
    cluster_corners = corners(close_corners_idx, :);
    merged_corner = mean(cluster_corners, 1);
    
    % 添加合并后的角点
    merged_corners = [merged_corners; merged_corner];
    merged(close_corners_idx) = true;
end
end 