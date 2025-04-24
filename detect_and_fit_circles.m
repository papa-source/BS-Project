function [centers, radii] = detect_and_fit_circles(edges)
% 圆孔识别与拟合
% 使用边缘检测结果进行圆孔检测和拟合

% 1. 预处理边缘图像
se = strel('disk', 1);
edges_cleaned = imclose(edges, se);
edges_cleaned = bwareaopen(edges_cleaned, 30);  % 增加面积阈值以去除更多噪声

% 2. 使用Hough变换进行圆检测
[centers, radii] = imfindcircles(edges_cleaned, [23 25], ... % 缩小半径范围以匹配实际圆孔大小
    'ObjectPolarity', 'dark', ...
    'Sensitivity', 0.92, ...  % 提高灵敏度
    'EdgeThreshold', 0.1);    % 降低边缘阈值以检测更多潜在圆

% 3. 如果没有检测到足够的圆，尝试使用连通区域分析
if size(centers, 1) < 4
    [L, num] = bwlabel(~edges_cleaned, 8);
    stats = regionprops(L, 'Area', 'Centroid', 'Perimeter', 'BoundingBox');
    
    for i = 1:num
        % 计算圆度
        circularity = 4 * pi * stats(i).Area / (stats(i).Perimeter^2);
        
        % 计算长宽比
        bbox = stats(i).BoundingBox;
        aspect_ratio = bbox(3) / bbox(4);
        
        % 严格的圆形筛选条件
        if circularity > 0.9 && ... % 提高圆度要求
            stats(i).Area > 100 && stats(i).Area < 500 && ... % 更严格的面积范围
            aspect_ratio > 0.95 && aspect_ratio < 1.05 % 更严格的长宽比限制
            
            % 获取当前区域的边缘点
            [y, x] = find(L == i);
            points = [x, y];
            
            % 使用最小二乘法拟合圆
            [center, radius] = fit_circle_to_points(points);
            
            if ~isempty(center) && radius > 10 && radius < 20  % 更严格的半径范围
                % 验证拟合结果
                quality = evaluate_circle_quality(edges, center, radius);
                if quality > 0.6  % 提高质量要求
                    centers = [centers; center];
                    radii = [radii; radius];
                end
            end
        end
    end
end

% 4. 移除重复的圆
if ~isempty(centers)
    valid_circles = true(size(centers, 1), 1);
    for i = 1:size(centers, 1)
        if ~valid_circles(i)
            continue;
        end
        % 计算当前圆与其他圆的距离
        distances = sqrt(sum((centers - centers(i,:)).^2, 2));
        overlaps = distances < 40;  % 增加距离阈值以更好地检测重复圆
        overlaps(i) = false;
        
        if any(overlaps)
            overlapping_idx = find(overlaps);
            qualities = zeros(length(overlapping_idx) + 1, 1);
            qualities(1) = evaluate_circle_quality(edges, centers(i,:), radii(i));
            for j = 1:length(overlapping_idx)
                qualities(j+1) = evaluate_circle_quality(edges, ...
                    centers(overlapping_idx(j),:), radii(overlapping_idx(j)));
            end
            [max_quality, best_idx] = max(qualities);
            % 只保留质量明显更好的圆
            if best_idx ~= 1 && max_quality > qualities(1) * 1.3  % 提高质量差异要求
                valid_circles(i) = false;
            else
                valid_circles(overlapping_idx) = false;
            end
        end
    end
    
    centers = centers(valid_circles, :);
    radii = radii(valid_circles);
end

% 5. 如果检测到的圆太多，只保留质量最好的6个
if size(centers, 1) > 6
    qualities = zeros(size(centers, 1), 1);
    for i = 1:size(centers, 1)
        qualities(i) = evaluate_circle_quality(edges, centers(i,:), radii(i));
    end
    [~, idx] = sort(qualities, 'descend');
    centers = centers(idx(1:6), :);
    radii = radii(idx(1:6));
end

% 6. 半径修正
radii = radii * 0.99;  % 轻微调整半径
end

function [center, radius] = fit_circle_to_points(points)
% 使用最小二乘法拟合圆
x = points(:,1);
y = points(:,2);

% 构建方程组
A = [-2*x, -2*y, ones(size(x))];
b = -(x.^2 + y.^2);

try
    % 求解最小二乘问题
    params = (A'*A)\(A'*b);
    
    % 提取圆心和半径
    center = [-params(1), -params(2)];
    radius = sqrt(params(1)^2 + params(2)^2 - params(3));
    
    % 验证结果
    if ~isreal(radius) || radius <= 0
        center = [];
        radius = [];
    end
catch
    center = [];
    radius = [];
end
end

function quality = evaluate_circle_quality(edges, center, radius)
% 评估圆的质量
theta = linspace(0, 2*pi, 120);  % 进一步增加采样点
x = round(center(1) + radius*cos(theta));
y = round(center(2) + radius*sin(theta));

% 确保坐标在图像范围内
valid_idx = x > 0 & x <= size(edges, 2) & y > 0 & y <= size(edges, 1);
x = x(valid_idx);
y = y(valid_idx);

if isempty(x) || isempty(y)
    quality = 0;
    return;
end

% 计算边缘点比例
idx = sub2ind(size(edges), y, x);
edge_points = sum(edges(idx));
quality = edge_points / length(idx);

% 检查内部区域
inner_radius = radius * 0.8;  % 扩大内部检查区域
x_inner = round(center(1) + inner_radius*cos(theta));
y_inner = round(center(2) + inner_radius*sin(theta));
valid_idx = x_inner > 0 & x_inner <= size(edges, 2) & ...
            y_inner > 0 & y_inner <= size(edges, 1);
x_inner = x_inner(valid_idx);
y_inner = y_inner(valid_idx);

if ~isempty(x_inner) && ~isempty(y_inner)
    idx_inner = sub2ind(size(edges), y_inner, x_inner);
    inner_points = sum(edges(idx_inner)) / length(idx_inner);
    if inner_points > 0.2  % 提高内部边缘点阈值
        quality = quality * 0.3;  % 加强对内部边缘点的惩罚
    end
end

% 考虑圆的完整性
if quality < 0.5  % 提高完整性要求
    quality = quality * 0.4;
end

% 考虑圆的对称性
outer_radius = radius * 1.1;
x_outer = round(center(1) + outer_radius*cos(theta));
y_outer = round(center(2) + outer_radius*sin(theta));
valid_idx = x_outer > 0 & x_outer <= size(edges, 2) & ...
            y_outer > 0 & y_outer <= size(edges, 1);
x_outer = x_outer(valid_idx);
y_outer = y_outer(valid_idx);

if ~isempty(x_outer) && ~isempty(y_outer)
    idx_outer = sub2ind(size(edges), y_outer, x_outer);
    outer_points = sum(edges(idx_outer)) / length(idx_outer);
    if outer_points > 0.1  % 如果外部有太多边缘点
        quality = quality * 0.5;
    end
end
end 