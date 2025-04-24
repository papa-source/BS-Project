function [centers, radii] = detect_and_fit_circles(edges)
% Բ��ʶ�������
% ʹ�ñ�Ե���������Բ�׼������

% 1. Ԥ�����Եͼ��
se = strel('disk', 1);
edges_cleaned = imclose(edges, se);
edges_cleaned = bwareaopen(edges_cleaned, 30);  % ���������ֵ��ȥ����������

% 2. ʹ��Hough�任����Բ���
[centers, radii] = imfindcircles(edges_cleaned, [23 25], ... % ��С�뾶��Χ��ƥ��ʵ��Բ�״�С
    'ObjectPolarity', 'dark', ...
    'Sensitivity', 0.92, ...  % ���������
    'EdgeThreshold', 0.1);    % ���ͱ�Ե��ֵ�Լ�����Ǳ��Բ

% 3. ���û�м�⵽�㹻��Բ������ʹ����ͨ�������
if size(centers, 1) < 4
    [L, num] = bwlabel(~edges_cleaned, 8);
    stats = regionprops(L, 'Area', 'Centroid', 'Perimeter', 'BoundingBox');
    
    for i = 1:num
        % ����Բ��
        circularity = 4 * pi * stats(i).Area / (stats(i).Perimeter^2);
        
        % ���㳤���
        bbox = stats(i).BoundingBox;
        aspect_ratio = bbox(3) / bbox(4);
        
        % �ϸ��Բ��ɸѡ����
        if circularity > 0.9 && ... % ���Բ��Ҫ��
            stats(i).Area > 100 && stats(i).Area < 500 && ... % ���ϸ�������Χ
            aspect_ratio > 0.95 && aspect_ratio < 1.05 % ���ϸ�ĳ��������
            
            % ��ȡ��ǰ����ı�Ե��
            [y, x] = find(L == i);
            points = [x, y];
            
            % ʹ����С���˷����Բ
            [center, radius] = fit_circle_to_points(points);
            
            if ~isempty(center) && radius > 10 && radius < 20  % ���ϸ�İ뾶��Χ
                % ��֤��Ͻ��
                quality = evaluate_circle_quality(edges, center, radius);
                if quality > 0.6  % �������Ҫ��
                    centers = [centers; center];
                    radii = [radii; radius];
                end
            end
        end
    end
end

% 4. �Ƴ��ظ���Բ
if ~isempty(centers)
    valid_circles = true(size(centers, 1), 1);
    for i = 1:size(centers, 1)
        if ~valid_circles(i)
            continue;
        end
        % ���㵱ǰԲ������Բ�ľ���
        distances = sqrt(sum((centers - centers(i,:)).^2, 2));
        overlaps = distances < 40;  % ���Ӿ�����ֵ�Ը��õؼ���ظ�Բ
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
            % ֻ�����������Ը��õ�Բ
            if best_idx ~= 1 && max_quality > qualities(1) * 1.3  % �����������Ҫ��
                valid_circles(i) = false;
            else
                valid_circles(overlapping_idx) = false;
            end
        end
    end
    
    centers = centers(valid_circles, :);
    radii = radii(valid_circles);
end

% 5. �����⵽��Բ̫�ֻ࣬����������õ�6��
if size(centers, 1) > 6
    qualities = zeros(size(centers, 1), 1);
    for i = 1:size(centers, 1)
        qualities(i) = evaluate_circle_quality(edges, centers(i,:), radii(i));
    end
    [~, idx] = sort(qualities, 'descend');
    centers = centers(idx(1:6), :);
    radii = radii(idx(1:6));
end

% 6. �뾶����
radii = radii * 0.99;  % ��΢�����뾶
end

function [center, radius] = fit_circle_to_points(points)
% ʹ����С���˷����Բ
x = points(:,1);
y = points(:,2);

% ����������
A = [-2*x, -2*y, ones(size(x))];
b = -(x.^2 + y.^2);

try
    % �����С��������
    params = (A'*A)\(A'*b);
    
    % ��ȡԲ�ĺͰ뾶
    center = [-params(1), -params(2)];
    radius = sqrt(params(1)^2 + params(2)^2 - params(3));
    
    % ��֤���
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
% ����Բ������
theta = linspace(0, 2*pi, 120);  % ��һ�����Ӳ�����
x = round(center(1) + radius*cos(theta));
y = round(center(2) + radius*sin(theta));

% ȷ��������ͼ��Χ��
valid_idx = x > 0 & x <= size(edges, 2) & y > 0 & y <= size(edges, 1);
x = x(valid_idx);
y = y(valid_idx);

if isempty(x) || isempty(y)
    quality = 0;
    return;
end

% �����Ե�����
idx = sub2ind(size(edges), y, x);
edge_points = sum(edges(idx));
quality = edge_points / length(idx);

% ����ڲ�����
inner_radius = radius * 0.8;  % �����ڲ��������
x_inner = round(center(1) + inner_radius*cos(theta));
y_inner = round(center(2) + inner_radius*sin(theta));
valid_idx = x_inner > 0 & x_inner <= size(edges, 2) & ...
            y_inner > 0 & y_inner <= size(edges, 1);
x_inner = x_inner(valid_idx);
y_inner = y_inner(valid_idx);

if ~isempty(x_inner) && ~isempty(y_inner)
    idx_inner = sub2ind(size(edges), y_inner, x_inner);
    inner_points = sum(edges(idx_inner)) / length(idx_inner);
    if inner_points > 0.2  % ����ڲ���Ե����ֵ
        quality = quality * 0.3;  % ��ǿ���ڲ���Ե��ĳͷ�
    end
end

% ����Բ��������
if quality < 0.5  % ���������Ҫ��
    quality = quality * 0.4;
end

% ����Բ�ĶԳ���
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
    if outer_points > 0.1  % ����ⲿ��̫���Ե��
        quality = quality * 0.5;
    end
end
end 