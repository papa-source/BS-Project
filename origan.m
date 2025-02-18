function hole_detection_system()
    % 板材冲孔质量检测系统
    % 包含图像预处理、边缘提取、圆孔识别与拟合、轮廓线识别与拟合
    
    % 读取图像
    [filename, pathname] = uigetfile({'*.jpg;*.png;*.bmp', '图像文件 (*.jpg, *.png, *.bmp)'});
    if filename == 0
        return;
    end
    img = imread(fullfile(pathname, filename));
    
    % 1. 图像预处理
    % 1.1 灰度化
    if size(img, 3) == 3
        gray_img = rgb2gray(img);
    else
        gray_img = img;
    end
    
    % 1.2 自适应中值滤波去噪
    filtered_img = adaptive_median_filter(gray_img);
    
    % 2. 图像分割
    % 2.1 迭代法求最佳阈值
    threshold = iterative_threshold(filtered_img);
    binary_img = filtered_img > threshold;
    % 2.2 ROI提取
roi = extract_roi(binary_img);

% 3. 边缘提取
% 使用Canny算子
edges = edge(roi, 'Canny');

% 4. 圆孔识别与拟合
[centers, radii] = detect_and_fit_circles(edges);
% 5. 轮廓线识别与拟合
[lines, corners] = detect_and_fit_contours(edges);

% 显示结果
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
    figure('Name', '板材冲孔质量检测结果');
    
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
    viscircles(centers, radii, 'EdgeColor', 'b');
end
title('圆孔检测结果');

% 轮廓线检测结果
subplot(2,2,4);
imshow(original_img);
hold on;
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1), xy(:,2), 'LineWidth', 2, 'Color', 'green');
end
if ~isempty(corners)
    plot(corners(:,1), corners(:,2), 'r*');
end
title('轮廓线检测结果');
end        
hasIPT = license('test', 'Image_Toolbox');
if ~hasIPT
    % 如果没有安装，会提示你需要安装 Image Processing Toolbox
    error('需要安装 Image Processing Toolbox');
end