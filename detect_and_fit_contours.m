function [lines, corners] = detect_and_fit_contours(edges)
% ������ʶ�������
% ֱ��ʹ�ñ�Ե�������������

% 1. �������������ڲ�Բ��
% ʹ����̬ѧ��������������
se = strel('disk', 2);
outer_mask = imclose(edges, se);
outer_edges = bwperim(outer_mask);

% ��ȡ�ڲ�Բ��
holes_mask = edges & ~outer_edges;
[L_holes, num_holes] = bwlabel(holes_mask, 8);
hole_stats = regionprops(L_holes, 'Centroid', 'Area', 'PixelList', 'Perimeter');

% 2. ����������
% ��ȡ�������ĵ㼯
[outer_y, outer_x] = find(outer_edges);
outer_points = [outer_x, outer_y];

% ��ʼ���߶νṹ
lines = struct('point1', {}, 'point2', {});

if ~isempty(outer_points)
    % ʹ��͹����ȡ��Ҫ�����㣬��ʹ�ø���ĵ��Ա�����״
    try
        k = boundary(outer_points(:,1), outer_points(:,2), 0.8); % ʹ�ý�С����������
        hull_points = outer_points(k,:);
        
        % �����������߶�
        for i = 1:length(k)-1
            lines(end+1).point1 = hull_points(i,:);
            lines(end).point2 = hull_points(i+1,:);
        end
        % �պ�������
        lines(end+1).point1 = hull_points(end,:);
        lines(end).point2 = hull_points(1,:);
    catch
        % ���boundaryʧ�ܣ�ʹ�ø��򵥵ķ���
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

% 3. ����Բ��
for i = 1:num_holes
    hole_pixels = hole_stats(i).PixelList;
    
    % ����Բ��
    circularity = 4 * pi * hole_stats(i).Area / (hole_stats(i).Perimeter^2);
    
    % ֻ����Բ�εĿ�
    if circularity > 0.6 && hole_stats(i).Area > 50
        % ��ȡԲ�ױ߽��
        [y, x] = find(L_holes == i);
        boundary_points = [x, y];
        
        % ʹ�ñ߽�㴴��Բ������
        try
            % ʹ�ñ߽���ٻ�ȡ����ı߽��
            k = boundary(boundary_points(:,1), boundary_points(:,2), 1);
            ordered_points = boundary_points(k,:);
            
            % ʹ�ù̶������ĵ�������Բ
            num_segments = 16;  % ���Ӷ���ʹԲ��ƽ��
            step = max(1, floor(length(k)/num_segments));
            
            for j = 1:num_segments
                idx1 = 1 + mod((j-1)*step, length(k));
                idx2 = 1 + mod(j*step, length(k));
                
                lines(end+1).point1 = ordered_points(idx1,:);
                lines(end).point2 = ordered_points(idx2,:);
            end
        catch
            % ����߽����ʧ�ܣ�ʹ�����ĵ�͹̶��뾶
            center = hole_stats(i).Centroid;
            radius = sqrt(hole_stats(i).Area / pi);
            
            % ����Բ������
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

% 4. ��ȡ�ǵ㣨��������������
corners = [];
if ~isempty(lines)
    % �����������߶ε�����
    num_outer_lines = length(k);
    outer_corners = get_corners(lines(1:num_outer_lines));
    corners = outer_corners;
end

end

function corners = get_corners(lines)
% ��ȡ���Ż��ǵ�
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
    
    % �ϲ�����Ľǵ�
    if ~isempty(corners)
        corners = merge_close_corners(corners);
    end
end
end

function merged_corners = merge_close_corners(corners)
% �ϲ�����Ľǵ�
dist_threshold = 10;
n = size(corners, 1);
merged = false(n, 1);
merged_corners = [];

for i = 1:n
    if merged(i)
        continue;
    end
    
    % �ҵ��뵱ǰ�ǵ�����������ǵ�
    distances = sqrt(sum((corners - corners(i,:)).^2, 2));
    close_corners_idx = distances < dist_threshold;
    
    % ������Щ�ǵ��ƽ��λ��
    merged_corner = mean(corners(close_corners_idx, :), 1);
    merged_corners = [merged_corners; merged_corner];
    merged(close_corners_idx) = true;
end
end

function corner = line_intersection(line1, line2)
% ������ֱ�߽���
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

% ��齻���Ƿ��������߶���
if is_point_on_line_segment([x,y], line1) && is_point_on_line_segment([x,y], line2)
    corner = [x, y];
else
    corner = [];
end
end

function on_segment = is_point_on_line_segment(point, line)
% �жϵ��Ƿ����߶���
x = point(1);
y = point(2);
x1 = line.point1(1);
y1 = line.point1(2);
x2 = line.point2(1);
y2 = line.point2(2);

% �����Ƿ����߶εķ�Χ��
on_segment = x >= min(x1,x2) - 1 && x <= max(x1,x2) + 1 && ...
             y >= min(y1,y2) - 1 && y <= max(y1,y2) + 1;
end 