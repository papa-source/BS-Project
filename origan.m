function hole_detection_system()
    % ��ĳ���������ϵͳ
    % ����ͼ��Ԥ������Ե��ȡ��Բ��ʶ������ϡ�������ʶ�������
    
    % ��ȡͼ��
    [filename, pathname] = uigetfile({'*.jpg;*.png;*.bmp', 'ͼ���ļ� (*.jpg, *.png, *.bmp)'});
    if filename == 0
        return;
    end
    img = imread(fullfile(pathname, filename));
    
    % 1. ͼ��Ԥ����
    % 1.1 �ҶȻ�
    if size(img, 3) == 3
        gray_img = rgb2gray(img);
    else
        gray_img = img;
    end
    
    % 1.2 ����Ӧ��ֵ�˲�ȥ��
    filtered_img = adaptive_median_filter(gray_img);
    
    % 2. ͼ��ָ�
    % 2.1 �������������ֵ
    threshold = iterative_threshold(filtered_img);
    binary_img = filtered_img > threshold;
    % 2.2 ROI��ȡ
roi = extract_roi(binary_img);

% 3. ��Ե��ȡ
% ʹ��Canny����
edges = edge(roi, 'Canny');

% 4. Բ��ʶ�������
[centers, radii] = detect_and_fit_circles(edges);
% 5. ������ʶ�������
[lines, corners] = detect_and_fit_contours(edges);

% ��ʾ���
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
    figure('Name', '��ĳ�����������');
    
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
    viscircles(centers, radii, 'EdgeColor', 'b');
end
title('Բ�׼����');

% �����߼����
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
title('�����߼����');
end        
hasIPT = license('test', 'Image_Toolbox');
if ~hasIPT
    % ���û�а�װ������ʾ����Ҫ��װ Image Processing Toolbox
    error('��Ҫ��װ Image Processing Toolbox');
end