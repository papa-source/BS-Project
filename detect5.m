function hole_detection_system()
% ��ĳ���������ϵͳ
% ����ͼ��Ԥ������Ե��ȡ��Բ��ʶ������ϡ�������ʶ��������ĸ���Ҫģ��

% ��ȡͼ��
[filename, pathname] = uigetfile({'*.jpg;*.png;*.bmp', 'ͼ���ļ� (*.jpg, *.png, *.bmp)'});
if filename == 0
    return;
end
img = imread(fullfile(pathname, filename));

%% 1. ͼ��Ԥ����ģ��
% 1.1 ͼ��ҶȻ�
if size(img, 3) == 3
    gray_img = rgb2gray(img);
else
    gray_img = img;
end

% 1.2 ����Ӧ��ֵ�˲�ȥ��
filtered_img = adaptive_median_filter(gray_img);

% 1.3 �Աȶ���ǿ
filtered_img = imadjust(filtered_img, stretchlim(filtered_img, [0.05 0.95]));

%% 2. ͼ��ָ���ROI��ȡģ��
% 2.1 ʹ�ö༶��ֵ�ָ�
T1 = graythresh(filtered_img);  % Otsu��������ȫ����ֵ
binary_img = imbinarize(filtered_img, T1);

% ȷ��Ŀ��Ϊ��ɫ��ֵΪ1��������Ϊ��ɫ��ֵΪ0��
if mean(binary_img(:)) < 0.5
    binary_img = ~binary_img;
end

% 2.2 ��̬ѧ������Ʊ�Ե
se1 = strel('disk', 3);
% �Ƚ��п�����ȥ�����
binary_img = imopen(binary_img, se1);
% �ٽ��б��������ӱ�Ե
binary_img = imclose(binary_img, se1);

% 2.3 �����ڲ�Բ��
% �����ͨ����
[L, num] = bwlabel(~binary_img, 8);
stats = regionprops(L, 'Area', 'Circularity');

% �ҳ�Բ������Բ�ף�
is_hole = false(num, 1);
for i = 1:num
    if stats(i).Area > 50 && stats(i).Area < 5000 && stats(i).Circularity > 0.8
        is_hole(i) = true;
    end
end

% �ؽ���ֵͼ�񣬱���Բ��
binary_img_with_holes = binary_img;
for i = 1:num
    if is_hole(i)
        binary_img_with_holes(L == i) = 0;  % ��Բ��������Ϊ��ɫ
    end
end

% 2.4 ROI��ȡ
% �Ƴ�С�������
binary_img = binary_img_with_holes;
binary_img = bwareaopen(binary_img, 50);  % ���������ֵ�Ա�������ϸ��
roi = extract_roi(binary_img);

%% 3. ��Ե��ȡģ��
% 3.1 Canny��Ե���
% ʹ�ö�ֵͼ����б�Ե��⣬�������˲����ͼ��
[~, thresh] = edge(binary_img, 'Canny');
edges = edge(binary_img, 'Canny', [0.1 0.3]);  % ʹ�ýϵ͵���ֵ�Լ������Ե

% 3.2 ����ڲ�Բ�ױ�Ե
% ʹ��Sobel������ǿԲ�α�Ե
[Gx, Gy] = imgradientxy(filtered_img, 'sobel');
grad_mag = imgradient(Gx, Gy);
circle_edges = grad_mag > mean(grad_mag(:)) + 2*std(grad_mag(:));
circle_edges = bwareaopen(circle_edges, 20);  % �Ƴ�С����

% 3.3 ��̬ѧ�����Ż���Ե
se1 = strel('line', 3, 0);    % ˮƽ��
se2 = strel('line', 3, 90);   % ��ֱ��
se3 = strel('line', 3, 45);   % 45����
se4 = strel('line', 3, -45);  % -45����
se5 = strel('disk', 2);       % Բ�νṹԪ�أ�������ǿԲ�ױ�Ե

% ���ӶϿ��ı�Ե
edges = imdilate(edges | circle_edges, [se1 se2 se3 se4 se5]);
edges = bwmorph(edges, 'bridge');  % ���ӶϿ��ı�Ե
edges = bwmorph(edges, 'clean');   % �Ƴ�������
edges = bwmorph(edges, 'spur');    % �Ƴ�ë��
edges = bwmorph(edges, 'thin', Inf);  % ϸ����Ե

% 3.4 ʹ���������Թ���С�ı�Ե��
edges = bwareaopen(edges, 30);  % �Ƴ�С��30���ص���ͨ����

% 3.5 ��ȡ�����������ڲ���Ե
% ��ȡ��ֵͼ��ı߽磬�����ڲ��׵ı߽�
boundary = bwperim(binary_img);
% �ϲ����б�Ե
edges = edges | boundary;

% ��ʾ�м������ڵ���
figure('Name', 'ͼ�����м���', 'NumberTitle', 'off');
subplot(3,2,1); imshow(gray_img); title('�Ҷ�ͼ��');
subplot(3,2,2); imshow(filtered_img); title('�˲��ͶԱȶ���ǿ');
subplot(3,2,3); imshow(binary_img); title('��ֵ�����');
subplot(3,2,4); imshow(roi); title('ROI��ȡ���');
subplot(3,2,5); 
imshow(edges);
title('��Ե�����������Բ�ף�');

% ��ʾֱ��ͼ
subplot(3,2,6); 
imhist(filtered_img);
title('ͼ��ֱ��ͼ');
hold on;
y_limits = ylim;
plot([T1*255 T1*255], [0 y_limits(2)], 'r-', 'LineWidth', 2);
legend('ֱ��ͼ', '��ֵ');
hold off;

%% 4. Բ��ʶ�������ģ��
% ʹ�ñ�Ե���������Բ���
[centers, radii] = imfindcircles(edges, [20 100], ...
    'ObjectPolarity', 'bright', ...
    'Sensitivity', 0.95, ...
    'EdgeThreshold', 0.1, ...
    'Method', 'PhaseCode');

% ���û�м�⵽Բ������ʹ�ö�ֵ����
if isempty(centers)
    holes_mask = ~binary_img_with_holes & binary_img;
    [centers, radii] = imfindcircles(holes_mask, [20 100], ...
        'ObjectPolarity', 'dark', ...
        'Sensitivity', 0.92, ...
        'EdgeThreshold', 0.1);
end

%% 5. ������ʶ�������ģ��
% ʹ�ù����Ķ�ֵͼ���Ե�����������
workpiece_edges = bwperim(binary_img_with_holes);

% ��ȡROI����
stats = regionprops(binary_img_with_holes, 'BoundingBox', 'Area');
if ~isempty(stats)
    % ѡ����������������Ϊ����
    [~, max_idx] = max([stats.Area]);
    bbox = stats(max_idx).BoundingBox;
    roi_x = round(bbox(1));
    roi_y = round(bbox(2));
    roi_w = round(bbox(3));
    roi_h = round(bbox(4));
    
    % ����ROI���룬����ӱ߾�
    margin = 10;  % ���10���صı߾�
    roi_mask = false(size(workpiece_edges));
    roi_y_min = max(1, roi_y - margin);
    roi_y_max = min(size(workpiece_edges,1), roi_y + roi_h + margin);
    roi_x_min = max(1, roi_x - margin);
    roi_x_max = min(size(workpiece_edges,2), roi_x + roi_w + margin);
    roi_mask(roi_y_min:roi_y_max, roi_x_min:roi_x_max) = true;
    
    % ֻ����ROI�����ڵı�Ե
    workpiece_edges = workpiece_edges & roi_mask;
    
    % �Ƴ�ͼ��߿�
    border_margin = 20;  % ���ñ߿�����Ŀ��
    workpiece_edges(1:border_margin,:) = 0;  % ����ϱ߿�
    workpiece_edges(end-border_margin:end,:) = 0;  % ����±߿�
    workpiece_edges(:,1:border_margin) = 0;  % �����߿�
    workpiece_edges(:,end-border_margin:end) = 0;  % ����ұ߿�
end

% ���������
[H_outer, theta_outer, rho_outer] = hough(workpiece_edges);
P_outer = houghpeaks(H_outer, 12, 'threshold', ceil(0.25*max(H_outer(:))));  % ���ӷ�ֵ��������������ֵ
lines = houghlines(workpiece_edges, theta_outer, rho_outer, P_outer, ...
    'FillGap', 25, 'MinLength', 25);  % ��������϶����С��С����

% �ϲ�������߶�
if ~isempty(lines)
    merged_lines = [];
    used = false(1, length(lines));
    
    for i = 1:length(lines)
        if used(i)
            continue;
        end
        
        current_line = lines(i);
        used(i) = true;
        
        % ���㵱ǰ�߶εĽǶ�
        dx = current_line.point2(1) - current_line.point1(1);
        dy = current_line.point2(2) - current_line.point1(2);
        angle1 = atan2(dy, dx);
        
        % Ѱ�ҿ��Ժϲ����߶�
        for j = i+1:length(lines)
            if used(j)
                continue;
            end
            
            % ������Ƚ��߶εĽǶ�
            dx = lines(j).point2(1) - lines(j).point1(1);
            dy = lines(j).point2(2) - lines(j).point1(2);
            angle2 = atan2(dy, dx);
            
            % �ſ�ϲ�����
            angle_diff = abs(mod(angle1 - angle2 + pi, pi));
            if angle_diff < 0.3 && ...  % ���ӽǶ��ݲ�
               (norm(current_line.point2 - lines(j).point1) < 40 || ...  % ���Ӿ����ݲ�
                norm(current_line.point1 - lines(j).point2) < 40)
                
                % �ϲ������߶�
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

% ��������պϴ���
if ~isempty(lines)
    % �ҵ������߶ζ˵�
    endpoints = [];
    for i = 1:length(lines)
        endpoints = [endpoints; lines(i).point1; lines(i).point2];
    end
    
    % Ѱ��δ�պϵĶ˵㣨�������˵�����Զ�ĵ㣩
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
        if min_dist > 30  % ���ñպ���ֵ
            unclosed(i) = true;
        end
    end
    
    % ����δ�պϵĶ˵�
    unclosed_points = endpoints(unclosed,:);
    if size(unclosed_points, 1) >= 2
        for i = 1:2:size(unclosed_points, 1)
            if i+1 <= size(unclosed_points, 1)
                % ������ԭʼ�߶νṹ��ͬ�����߶�
                new_line = struct('point1', unclosed_points(i,:), ...
                                'point2', unclosed_points(i+1,:), ...
                                'theta', 0, ...
                                'rho', 0);
                lines = [lines, new_line];
            end
        end
    end
end

% ��ȡ�ǵ�
corners = [];
if length(lines) >= 2
    for i = 1:length(lines)-1
        for j = i+1:length(lines)
            corner = line_intersection(lines(i), lines(j));
            if ~isempty(corner)
                % ���ǵ��Ƿ���ROI������
                if corner(1) >= roi_x && corner(1) <= roi_x+roi_w && ...
                   corner(2) >= roi_y && corner(2) <= roi_y+roi_h
                    corners = [corners; corner];
                end
            end
        end
    end
    
    % �ϲ�����Ľǵ�
    if ~isempty(corners)
        corners = merge_close_corners(corners);
    end
end

%% 6. ��ʾ���
display_results(img, edges, centers, radii, lines, corners);
end

function filtered = adaptive_median_filter(img)
% ����Ӧ��ֵ�˲�ʵ��
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
% ��ȡͼ�񴰿�
half = floor(window_size/2);
[m, n] = size(img);
row_min = max(1, i-half);
row_max = min(m, i+half);
col_min = max(1, j-half);
col_max = min(n, j+half);
window = img(row_min:row_max, col_min:col_max);
end

function threshold = iterative_threshold(img)
% �������������ֵ
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
% ROI��ȡ
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
% Բ��ʶ�������
[centers, radii] = imfindcircles(edges, [20 100], 'Sensitivity', 0.85);
if isempty(centers)
    centers = [];
    radii = [];
end
end

function [lines, corners] = detect_and_fit_contours(edges)
% ������ʶ�������
[H, theta, rho] = hough(edges);
P = houghpeaks(H, 10, 'threshold', ceil(0.3*max(H(:))));
lines = houghlines(edges, theta, rho, P, 'FillGap', 5, 'MinLength', 7);

% ��ȡ�ǵ�
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
corner = [x, y];
end

function display_results(original_img, edges, centers, radii, lines, corners)
% ��ʾ�����
figure('Name', '��ĳ�����������', 'NumberTitle', 'off');

% ԭͼ
subplot(2,2,1);
imshow(original_img);
title('ԭʼͼ��');

% ��Ե�����
subplot(2,2,2);
imshow(edges);
title('��Ե�����');

% Բ�׼����
subplot(2,2,3);
imshow(original_img);
hold on;
if ~isempty(centers)
    viscircles(centers, radii, 'EdgeColor', 'b', 'LineWidth', 1.5);
    % ��עԲ�ĺͰ뾶
    for i = 1:size(centers, 1)
        plot(centers(i,1), centers(i,2), 'b+', 'MarkerSize', 10, 'LineWidth', 2);
        text(centers(i,1)+10, centers(i,2)+10, ...
             sprintf('R=%.1f', radii(i)), ...
             'Color', 'blue', 'FontSize', 8);
    end
end
title('Բ�׼����');
hold off;

% �����߼����
subplot(2,2,4);
imshow(original_img);
hold on;
% ����������
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1), xy(:,2), 'LineWidth', 2, 'Color', 'green');
end
% ���ƽǵ�
if ~isempty(corners)
    plot(corners(:,1), corners(:,2), 'r*', 'MarkerSize', 10);
end
title('�����߼����');
hold off;
end

function merged_corners = merge_close_corners(corners)
% �ϲ�����Ľǵ�
if isempty(corners)
    merged_corners = corners;
    return;
end

% ���þ�����ֵ
dist_threshold = 10;

% ��ʼ���������
n = size(corners, 1);
merged = false(n, 1);
merged_corners = [];

% �������нǵ�
for i = 1:n
    if merged(i)
        continue;
    end
    
    % ���㵱ǰ�ǵ��������ǵ�ľ���
    distances = sqrt(sum((corners - corners(i,:)).^2, 2));
    close_corners_idx = distances < dist_threshold;
    
    % ��������ǵ��ƽ��λ��
    cluster_corners = corners(close_corners_idx, :);
    merged_corner = mean(cluster_corners, 1);
    
    % ��Ӻϲ���Ľǵ�
    merged_corners = [merged_corners; merged_corner];
    merged(close_corners_idx) = true;
end
end 