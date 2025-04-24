function [lines, corners] = detect_and_fit_contours(edges)
% 轮廓线识别与拟合
% 直接使用边缘检测结果构建轮廓

% 1. 分离外轮廓和内部圆孔
% 使用形态学操作分离外轮廓
se = strel('disk', 2);
outer_mask = imclose(edges, se);
outer_edges = bwperim(outer_mask);

% 获取内部圆孔
holes_mask = edges & ~outer_edges;
[L_holes, num_holes] = bwlabel(holes_mask, 8);
hole_stats = regionprops(L_holes, 'Centroid', 'Area', 'PixelList', 'Perimeter');

% 2. 处理外轮廓
% 获取外轮廓的点集
[outer_y, outer_x] = find(outer_edges);
outer_points = [outer_x, outer_y];

% 初始化线段结构
lines = struct('point1', {}, 'point2', {});

if ~isempty(outer_points)
    % 使用凸包获取主要轮廓点，但使用更多的点以保持形状
    try
        k = boundary(outer_points(:,1), outer_points(:,2), 0.8); % 使用较小的收缩因子
        hull_points = outer_points(k,:);
        
        % 创建外轮廓线段
        for i = 1:length(k)-1
            lines(end+1).point1 = hull_points(i,:);
            lines(end).point2 = hull_points(i+1,:);
        end
        % 闭合外轮廓
        lines(end+1).point1 = hull_points(end,:);
        lines(end).point2 = hull_points(1,:);
    catch
        % 如果boundary失败，使用更简单的方法
        k = convhull(outer_points(:,1), outer_points(:,2));
        hull_points = outer_points(k,:);
        
        for i = 1:length(k)-1
            lines(end+1).point1 = hull_points(i,:);
            lines(end).point2 = hull_points(i+1,:);
        end
        lines(end+1).point1 = hull_points(end,:);
        lines(end).point2 = hull_points(1,:);
    end
end

% 3. 处理圆孔
for i = 1:num_holes
    hole_pixels = hole_stats(i).PixelList;
    
    % 计算圆度
    circularity = 4 * pi * hole_stats(i).Area / (hole_stats(i).Perimeter^2);
    
    % 只处理圆形的孔
    if circularity > 0.6 && hole_stats(i).Area > 50
        % 获取圆孔边界点
        [y, x] = find(L_holes == i);
        boundary_points = [x, y];
        
        % 使用边界点创建圆形轮廓
        try
            % 使用边界跟踪获取有序的边界点
            k = boundary(boundary_points(:,1), boundary_points(:,2), 1);
            ordered_points = boundary_points(k,:);
            
            % 使用固定数量的点来近似圆
            num_segments = 16;  % 增加段数使圆更平滑
            step = max(1, floor(length(k)/num_segments));
            
            for j = 1:num_segments
                idx1 = 1 + mod((j-1)*step, length(k));
                idx2 = 1 + mod(j*step, length(k));
                
                lines(end+1).point1 = ordered_points(idx1,:);
                lines(end).point2 = ordered_points(idx2,:);
            end
        catch
            % 如果边界跟踪失败，使用中心点和固定半径
            center = hole_stats(i).Centroid;
            radius = sqrt(hole_stats(i).Area / pi);
            
            % 创建圆形轮廓
            theta = linspace(0, 2*pi, num_segments+1);
            for j = 1:num_segments
                lines(end+1).point1 = [center(1) + radius*cos(theta(j)), ...
                                     center(2) + radius*sin(theta(j))];
                lines(end).point2 = [center(1) + radius*cos(theta(j+1)), ...
                                     center(2) + radius*sin(theta(j+1))];
            end
        end
    end
end

% 4. 获取角点（仅处理外轮廓）
corners = [];
if ~isempty(lines)
    % 计算外轮廓线段的数量
    num_outer_lines = length(k);
    outer_corners = get_corners(lines(1:num_outer_lines));
    corners = outer_corners;
end

end

function corners = get_corners(lines)
% 获取并优化角点
corners = [];
if length(lines) >= 2
    for i = 1:length(lines)-1
        for j = i+1:length(lines)
            corner = line_intersection(lines(i), lines(j));
            if ~isempty(corner)
                corners = [corners; corner];
            end
        end
    end
    
    % 合并相近的角点
    if ~isempty(corners)
        corners = merge_close_corners(corners);
    end
end
end

function merged_corners = merge_close_corners(corners)
% 合并相近的角点
dist_threshold = 10;
n = size(corners, 1);
merged = false(n, 1);
merged_corners = [];

for i = 1:n
    if merged(i)
        continue;
    end
    
    % 找到与当前角点相近的其他角点
    distances = sqrt(sum((corners - corners(i,:)).^2, 2));
    close_corners_idx = distances < dist_threshold;
    
    % 计算这些角点的平均位置
    merged_corner = mean(corners(close_corners_idx, :), 1);
    merged_corners = [merged_corners; merged_corner];
    merged(close_corners_idx) = true;
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

% 检查交点是否在两条线段上
if is_point_on_line_segment([x,y], line1) && is_point_on_line_segment([x,y], line2)
    corner = [x, y];
else
    corner = [];
end
end

function on_segment = is_point_on_line_segment(point, line)
% 判断点是否在线段上
x = point(1);
y = point(2);
x1 = line.point1(1);
y1 = line.point1(2);
x2 = line.point2(1);
y2 = line.point2(2);

% 检查点是否在线段的范围内
on_segment = x >= min(x1,x2) - 1 && x <= max(x1,x2) + 1 && ...
             y >= min(y1,y2) - 1 && y <= max(y1,y2) + 1;
end 