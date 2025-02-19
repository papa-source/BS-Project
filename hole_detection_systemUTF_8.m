function hole_detection_system()
% 鏉挎潗鍐插瓟璐ㄩ噺妫�娴嬬郴缁?
% 鍖呭惈鍥惧儚棰勫鐞嗐�佽竟缂樻彁鍙栥�佸渾瀛旇瘑鍒笌鎷熷悎銆佽疆寤撶嚎璇嗗埆涓庢嫙鍚堝洓涓富瑕佹ā鍧?

% 璇诲彇鍥惧儚
[filename, pathname] = uigetfile({'*.jpg;*.png;*.bmp', '鍥惧儚鏂囦欢 (*.jpg, *.png, *.bmp)'});
if filename == 0
    return;
end
img = imread(fullfile(pathname, filename));

%% 1. 鍥惧儚棰勫鐞嗘ā鍧?
% 1.1 鍥惧儚鐏板害鍖栫洏
if size(img, 3) == 3
    gray_img = rgb2gray(img);
else
    gray_img = img;
endC:\Users\YTONY\Desktop\BS

% 1.2 鑷�傚簲涓�兼护娉㈠幓鍣?
filtered_img = adaptive_median_filter(gray_img);

% 1.3 瀵规瘮搴﹀寮?
filtered_img = imadjust(filtered_img, stretchlim(filtered_img, [0.05 0.95]));

%% 2. 鍥惧儚鍒嗗壊涓嶳OI鎻愬彇妯″潡
% 2.1 浣跨敤澶氱骇闃堝�煎垎鍓?
T1 = graythresh(filtered_img);  % Otsu鏂规硶璁＄畻鍏ㄥ眬闃堝�?
binary_img = imbinarize(filtered_img, T1);

% 纭繚鐩爣涓虹櫧鑹诧紙鍊间负1锛夛紝鑳屾櫙涓洪粦鑹诧紙鍊间负0锛?
if mean(binary_img(:)) < 0.5
    binary_img = ~binary_img;
end

% 2.2 褰㈡�佸澶勭悊鏀瑰杽杈圭紭
se1 = strel('disk', 2);  % 鍑忓皬缁撴瀯鍏冪礌澶у皬
% 鍏堣繘琛屽紑杩愮畻鍘婚櫎鍣偣
binary_img = imopen(binary_img, se1);
% 鍐嶈繘琛岄棴杩愮畻杩炴帴杈圭紭
binary_img = imclose(binary_img, se1);

% 2.3 淇濈暀鍐呴儴鍦嗗瓟
% 鏍囪杩為�氬尯鍩?
[L, num] = bwlabel(~binary_img, 8);
stats = regionprops(L, 'Area', 'Circularity');

% 鎵惧嚭鍦嗗舰鍖哄煙锛堝渾瀛旓級
is_hole = false(num, 1);
for i = 1:num
    if stats(i).Area > 30 && stats(i).Area < 1000 && stats(i).Circularity > 0.85  % 璋冩暣闈㈢Н鑼冨洿鍜屽渾搴﹂槇鍊?
        is_hole(i) = true;
    end
end

% 閲嶅缓浜屽�煎浘鍍忥紝淇濈暀鍦嗗瓟
binary_img_with_holes = binary_img;
for i = 1:num
    if is_hole(i)
        binary_img_with_holes(L == i) = 0;  % 灏嗗渾瀛斿尯鍩熻涓洪粦鑹?
    end
end

% 2.4 ROI鎻愬彇
% 绉婚櫎灏忛潰绉尯鍩?
binary_img = binary_img_with_holes;
binary_img = bwareaopen(binary_img, 50);  % 闄嶄綆闈㈢Н闃堝�间互淇濈暀鏇村缁嗚妭
roi = extract_roi(binary_img);

%% 3. 杈圭紭鎻愬彇妯″潡
% 3.1 浣跨敤浜屽�煎浘鍍忕洿鎺ユ彁鍙栬竟缂?
edges = bwperim(~binary_img_with_holes);  % 浣跨敤甯﹀瓟鐨勪簩鍊煎浘鍍忕洿鎺ユ彁鍙栬竟缂?

% 3.2 浼樺寲杈圭紭
se = strel('disk', 2);  % 澧炲姞缁撴瀯鍏冪礌澶у皬
edges = imdilate(edges, se);  % 杞诲井鑶ㄨ儉浠ュ寮鸿竟缂?
edges = bwmorph(edges, 'thin', Inf);  % 缁嗗寲杈圭紭
edges = bwareaopen(edges, 30);  % 澧炲姞闈㈢Н闃堝�硷紝鍘婚櫎鏇村鍣０

% 鏄剧ず涓棿缁撴灉鐢ㄤ簬璋冭瘯
figure('Name', '鍥惧儚澶勭悊涓棿缁撴灉', 'NumberTitle', 'off');
subplot(3,2,1); imshow(gray_img); title('鐏板害鍥惧儚');
subplot(3,2,2); imshow(filtered_img); title('婊ゆ尝鍜屽姣斿害澧炲己');
subplot(3,2,3); imshow(binary_img); title('浜屽�煎寲缁撴灉');
subplot(3,2,4); imshow(roi); title('ROI鎻愬彇缁撴灉');
subplot(3,2,5); 
imshow(edges);
title('杈圭紭妫�娴嬬粨鏋滐紙鍖呭惈鍦嗗瓟锛?');

% 鏄剧ず鐩存柟鍥?
subplot(3,2,6); 
imhist(filtered_img);
title('鍥惧儚鐩存柟鍥?');
hold on;
y_limits = ylim;
plot([T1*255 T1*255], [0 y_limits(2)], 'r-', 'LineWidth', 2);
legend('鐩存柟鍥?', '闃堝�?');
hold off;

%% 4. 鍦嗗瓟璇嗗埆涓庢嫙鍚堟ā鍧?
% 浣跨敤杈圭紭妫�娴嬬粨鏋滆繘琛屽渾妫�娴?
[centers, radii] = imfindcircles(edges, [10 30], ...  % 璋冩暣鍗婂緞鑼冨洿涓?10-30鍍忕礌
    'ObjectPolarity', 'bright', ...
    'Sensitivity', 0.85, ...  % 闄嶄綆鏁忔劅搴?
    'EdgeThreshold', 0.1, ...
    'Method', 'PhaseCode');

% 濡傛灉娌℃湁妫�娴嬪埌鍦嗭紝灏濊瘯浣跨敤浜屽�兼帺鐮?
if isempty(centers)
    holes_mask = ~binary_img_with_holes & binary_img;
    [centers, radii] = imfindcircles(holes_mask, [10 30], ...  % 杩欓噷涔熻璋冩暣鍗婂緞鑼冨洿
        'Sensitivity', 0.85, ...
        'EdgeThreshold', 0.1);
end

%% 5. 杞粨绾胯瘑鍒笌鎷熷悎妯″潡
% 浣跨敤宸ヤ欢鐨勪簩鍊煎浘鍍忚竟缂樿繘琛岃疆寤撴娴?
workpiece_edges = bwperim(binary_img_with_holes);

% 鑾峰彇ROI鍖哄煙
stats = regionprops(binary_img_with_holes, 'BoundingBox', 'Area');
if ~isempty(stats)
    % 閫夋嫨鏈�澶ч潰绉殑鍖哄煙浣滀负宸ヤ欢
    [~, max_idx] = max([stats.Area]);
    bbox = stats(max_idx).BoundingBox;
    roi_x = round(bbox(1));
    roi_y = round(bbox(2));
    roi_w = round(bbox(3));
    roi_h = round(bbox(4));
    
    % 鍒涘缓ROI鎺╃爜锛屽苟娣诲姞杈硅窛
    margin = 10;  % 娣诲姞10鍍忕礌鐨勮竟璺?
    roi_mask = false(size(workpiece_edges));
    roi_y_min = max(1, roi_y - margin);
    roi_y_max = min(size(workpiece_edges,1), roi_y + roi_h + margin);
    roi_x_min = max(1, roi_x - margin);
    roi_x_max = min(size(workpiece_edges,2), roi_x + roi_w + margin);
    roi_mask(roi_y_min:roi_y_max, roi_x_min:roi_x_max) = true;
    
    % 鍙繚鐣橰OI鍖哄煙鍐呯殑杈圭紭
    workpiece_edges = workpiece_edges & roi_mask;
    
    % 绉婚櫎鍥惧儚杈规
    border_margin = 20;  % 璁剧疆杈规鍖哄煙鐨勫搴?
    workpiece_edges(1:border_margin,:) = 0;  % 娓呴櫎涓婅竟妗?
    workpiece_edges(end-border_margin:end,:) = 0;  % 娓呴櫎涓嬭竟妗?
    workpiece_edges(:,1:border_margin) = 0;  % 娓呴櫎宸﹁竟妗?
    workpiece_edges(:,end-border_margin:end) = 0;  % 娓呴櫎鍙宠竟妗?
end

% 澶栬疆寤撴娴?
[H_outer, theta_outer, rho_outer] = hough(workpiece_edges);
P_outer = houghpeaks(H_outer, 12, 'threshold', ceil(0.25*max(H_outer(:))));  % 澧炲姞宄板�肩偣鏁伴噺锛岄檷浣庨槇鍊?
lines = houghlines(workpiece_edges, theta_outer, rho_outer, P_outer, ...
    'FillGap', 25, 'MinLength', 25);  % 澧炲姞濉厖闂撮殭锛屽噺灏忔渶灏忛暱搴?

% 鍚堝苟鐩歌繎鐨勭嚎娈?
if ~isempty(lines)
    merged_lines = [];
    used = false(1, length(lines));
    
    for i = 1:length(lines)
        if used(i)
            continue;
        end
        
        current_line = lines(i);
        used(i) = true;
        
        % 璁＄畻褰撳墠绾挎鐨勮搴?
        dx = current_line.point2(1) - current_line.point1(1);
        dy = current_line.point2(2) - current_line.point1(2);
        angle1 = atan2(dy, dx);
        
        % 瀵绘壘鍙互鍚堝苟鐨勭嚎娈?
        for j = i+1:length(lines)
            if used(j)
                continue;
            end
            
            % 璁＄畻寰呮瘮杈冪嚎娈电殑瑙掑害
            dx = lines(j).point2(1) - lines(j).point1(1);
            dy = lines(j).point2(2) - lines(j).point1(2);
            angle2 = atan2(dy, dx);
            
            % 鏀惧鍚堝苟鏉′欢
            angle_diff = abs(mod(angle1 - angle2 + pi, pi));
            if angle_diff < 0.3 && ...  % 澧炲姞瑙掑害瀹瑰樊
               (norm(current_line.point2 - lines(j).point1) < 40 || ...  % 澧炲姞璺濈瀹瑰樊
                norm(current_line.point1 - lines(j).point2) < 40)
                
                % 鍚堝苟涓ゆ潯绾挎
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

% 娣诲姞杞粨闂悎澶勭悊
if ~isempty(lines)
    % 鎵惧埌鎵�鏈夌嚎娈电鐐?
    endpoints = [];
    for i = 1:length(lines)
        endpoints = [endpoints; lines(i).point1; lines(i).point2];
    end
    
    % 瀵绘壘鏈棴鍚堢殑绔偣锛堜笌鍏朵粬绔偣璺濈杈冭繙鐨勭偣锛?
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
        if min_dist > 30  % 璁剧疆闂悎闃堝�?
            unclosed(i) = true;
        end
    end
    
    % 杩炴帴鏈棴鍚堢殑绔偣
    unclosed_points = endpoints(unclosed,:);
    if size(unclosed_points, 1) >= 2
        for i = 1:2:size(unclosed_points, 1)
            if i+1 <= size(unclosed_points, 1)
                % 鍒涘缓涓庡師濮嬬嚎娈电粨鏋勭浉鍚岀殑鏂扮嚎娈?
                new_line = struct('point1', unclosed_points(i,:), ...
                                'point2', unclosed_points(i+1,:), ...
                                'theta', 0, ...
                                'rho', 0);
                lines = [lines, new_line];
            end
        end
    end
end

% 鑾峰彇瑙掔偣
corners = [];
if length(lines) >= 2
    for i = 1:length(lines)-1
        for j = i+1:length(lines)
            corner = line_intersection(lines(i), lines(j));
            if ~isempty(corner)
                % 妫�鏌ヨ鐐规槸鍚﹀湪ROI鍖哄煙鍐?
                if corner(1) >= roi_x && corner(1) <= roi_x+roi_w && ...
                   corner(2) >= roi_y && corner(2) <= roi_y+roi_h
                    corners = [corners; corner];
                end
            end
        end
    end
    
    % 鍚堝苟鐩歌繎鐨勮鐐?
    if ~isempty(corners)
        corners = merge_close_corners(corners);
    end
end

%% 6. 鏄剧ず缁撴灉
display_results(img, edges, centers, radii, lines, corners);
end

function filtered = adaptive_median_filter(img)
% 鑷�傚簲涓�兼护娉㈠疄鐜?
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
% 鑾峰彇鍥惧儚绐楀彛
half = floor(window_size/2);
[m, n] = size(img);
row_min = max(1, i-half);
row_max = min(m, i+half);
col_min = max(1, j-half);
col_max = min(n, j+half);
window = img(row_min:row_max, col_min:col_max);
end

function threshold = iterative_threshold(img)
% 杩唬娉曟眰鏈�浣抽槇鍊?
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
% ROI鎻愬彇
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
% 鍦嗗瓟璇嗗埆涓庢嫙鍚?
[centers, radii] = imfindcircles(edges, [20 100], 'Sensitivity', 0.85);
if isempty(centers)
    centers = [];
    radii = [];
end
end

function [lines, corners] = detect_and_fit_contours(edges)
% 杞粨绾胯瘑鍒笌鎷熷悎
[H, theta, rho] = hough(edges);
P = houghpeaks(H, 10, 'threshold', ceil(0.3*max(H(:))));
lines = houghlines(edges, theta, rho, P, 'FillGap', 5, 'MinLength', 7);

% 鑾峰彇瑙掔偣
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
% 璁＄畻涓ょ洿绾夸氦鐐?
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
% 鏄剧ず妫�娴嬬粨鏋?
figure('Name', '鏉挎潗鍐插瓟璐ㄩ噺妫�娴嬬粨鏋?', 'NumberTitle', 'off');

% 鍘熷浘
subplot(2,2,1);
imshow(original_img);
title('鍘熷鍥惧儚');

% 杈圭紭妫�娴嬬粨鏋?
subplot(2,2,2);
imshow(edges);
title('杈圭紭妫�娴嬬粨鏋?');

% 鍦嗗瓟妫�娴嬬粨鏋?
subplot(2,2,3);
imshow(original_img);
hold on;
if ~isempty(centers)
    viscircles(centers, radii, 'EdgeColor', 'b', 'LineWidth', 1.5);
    % 鏍囨敞鍦嗗績鍜屽崐寰?
    for i = 1:size(centers, 1)
        plot(centers(i,1), centers(i,2), 'b+', 'MarkerSize', 10, 'LineWidth', 2);
        text(centers(i,1)+10, centers(i,2)+10, ...
             sprintf('R=%.1f', radii(i)), ...
             'Color', 'blue', 'FontSize', 8);
    end
end
title('鍦嗗瓟妫�娴嬬粨鏋?');
hold off;

% 杞粨绾挎娴嬬粨鏋?
subplot(2,2,4);
imshow(original_img);
hold on;
% 缁樺埗杞粨绾?
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1), xy(:,2), 'LineWidth', 2, 'Color', 'green');
end
% 缁樺埗瑙掔偣
if ~isempty(corners)
    plot(corners(:,1), corners(:,2), 'r*', 'MarkerSize', 10);
end
title('杞粨绾挎娴嬬粨鏋?');
hold off;
end

function merged_corners = merge_close_corners(corners)
% 鍚堝苟鐩歌繎鐨勮鐐?
if isempty(corners)
    merged_corners = corners;
    return;
end

% 璁剧疆璺濈闃堝�?
dist_threshold = 10;

% 鍒濆鍖栨爣璁版暟缁?
n = size(corners, 1);
merged = false(n, 1);
merged_corners = [];

% 閬嶅巻鎵�鏈夎鐐?
for i = 1:n
    if merged(i)
        continue;
    end
    
    % 璁＄畻褰撳墠瑙掔偣涓庡叾浠栬鐐圭殑璺濈
    distances = sqrt(sum((corners - corners(i,:)).^2, 2));
    close_corners_idx = distances < dist_threshold;
    
    % 璁＄畻鐩歌繎瑙掔偣鐨勫钩鍧囦綅缃?
    cluster_corners = corners(close_corners_idx, :);
    merged_corner = mean(cluster_corners, 1);
    
    % 娣诲姞鍚堝苟鍚庣殑瑙掔偣
    merged_corners = [merged_corners; merged_corner];
    merged(close_corners_idx) = true;
end
end 