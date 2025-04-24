function results = compare_detected_vs_design(detected_centers, detected_radii, design_centers, design_radii, pixel_scale)
% COMPARE_DETECTED_VS_DESIGN 比较检测到的圆孔与设计圆孔之间的差异
%
% 输入参数:
%   detected_centers - 检测到的圆孔中心坐标 [x,y] 的Nx2矩阵
%   detected_radii - 检测到的圆孔半径的Nx1向量
%   design_centers - DXF文件中的设计圆孔中心坐标 [x,y] 的Mx2矩阵
%   design_radii - DXF文件中的设计圆孔半径的Mx1向量
%   pixel_scale - 像素到设计单位的比例因子(mm/pixel)，默认为1
%
% 输出:
%   results - 包含比较结果的结构体，具有以下字段:
%     .matched_pairs - 匹配的检测圆孔和设计圆孔的索引对
%     .center_errors - 每对匹配圆孔的中心位置误差(设计单位)
%     .radius_errors - 每对匹配圆孔的半径误差(设计单位)
%     .missing_design - 未被检测到的设计圆孔索引
%     .extra_detected - 检测到但不在设计中的圆孔索引
%     .summary - 包含整体统计信息的子结构体

% 设置默认缩放比例
if nargin < 5
    pixel_scale = 1; % 默认像素到设计单位的比例
    disp('未提供像素比例，使用默认值1.0');
end

% 初始化输出结构体
results = struct('matched_pairs', [], 'center_errors', [], ...
    'radius_errors', [], 'missing_design', [], 'extra_detected', [], ...
    'summary', struct());

% 检查输入是否为空
if isempty(detected_centers) || isempty(design_centers)
    disp('检测到的圆孔或设计圆孔为空，无法进行比较');
    return;
end

% 将检测到的坐标转换为设计单位
detected_centers_scaled = detected_centers * pixel_scale;
detected_radii_scaled = detected_radii * pixel_scale;

% 初始化匹配矩阵和距离矩阵
n_detected = size(detected_centers_scaled, 1);
n_design = size(design_centers, 1);
dist_matrix = zeros(n_detected, n_design);

% 计算每个检测圆与每个设计圆之间的欧几里得距离
for i = 1:n_detected
    for j = 1:n_design
        dist_matrix(i,j) = norm(detected_centers_scaled(i,:) - design_centers(j,:));
    end
end

% 匹配圆孔 - 使用匈牙利算法进行最优匹配
% 但我们需要设置一个最大距离阈值，超过此距离的圆不应该被匹配
max_dist_threshold = mean([detected_radii_scaled; design_radii]) * 0.5; % 半径平均值的50%作为阈值

% 将超过阈值的距离设置为一个大数，确保不会匹配
dist_matrix_thresholded = dist_matrix;
dist_matrix_thresholded(dist_matrix > max_dist_threshold) = 1e6;

% 使用匈牙利算法找到最优匹配
if exist('assignDetectionsToTracks', 'file') == 2
    % 如果有Computer Vision Toolbox，使用其函数
    [assignments, unassigned_detections, unassigned_designs] = ...
        assignDetectionsToTracks(dist_matrix_thresholded');
    
    matched_pairs = [assignments(:,2), assignments(:,1)]; % [检测索引, 设计索引]
    extra_detected = unassigned_detections;
    missing_design = unassigned_designs;
else
    % 简单的贪婪匹配方法
    matched_pairs = [];
    matched_detected = false(n_detected, 1);
    matched_design = false(n_design, 1);
    
    % 对距离矩阵进行排序
    [sorted_dists, idx] = sort(dist_matrix_thresholded(:));
    [det_idx, des_idx] = ind2sub(size(dist_matrix_thresholded), idx);
    
    for i = 1:length(sorted_dists)
        if sorted_dists(i) >= 1e6
            break; % 超过阈值，不再匹配
        end
        
        d_idx = det_idx(i);
        g_idx = des_idx(i);
        
        if ~matched_detected(d_idx) && ~matched_design(g_idx)
            matched_pairs = [matched_pairs; d_idx, g_idx];
            matched_detected(d_idx) = true;
            matched_design(g_idx) = true;
        end
    end
    
    % 找出未匹配的检测和设计
    extra_detected = find(~matched_detected);
    missing_design = find(~matched_design);
end

% 计算每对匹配圆孔的误差
center_errors = zeros(size(matched_pairs, 1), 1);
radius_errors = zeros(size(matched_pairs, 1), 1);

for i = 1:size(matched_pairs, 1)
    det_idx = matched_pairs(i, 1);
    des_idx = matched_pairs(i, 2);
    
    % 计算中心位置误差
    center_errors(i) = norm(detected_centers_scaled(det_idx,:) - design_centers(des_idx,:));
    
    % 计算半径误差
    radius_errors(i) = detected_radii_scaled(det_idx) - design_radii(des_idx);
end

% 统计摘要
summary = struct();
summary.num_matched = size(matched_pairs, 1);
summary.num_missing = length(missing_design);
summary.num_extra = length(extra_detected);
summary.match_rate = summary.num_matched / n_design;
summary.mean_center_error = mean(center_errors);
summary.std_center_error = std(center_errors);
summary.mean_radius_error = mean(radius_errors);
summary.std_radius_error = std(radius_errors);

% 填充结果结构体
results.matched_pairs = matched_pairs;
results.center_errors = center_errors;
results.radius_errors = radius_errors;
results.missing_design = missing_design;
results.extra_detected = extra_detected;
results.summary = summary;

% 可视化比较结果
visualize_comparison(detected_centers_scaled, detected_radii_scaled, ...
    design_centers, design_radii, matched_pairs, missing_design, extra_detected);

end

function visualize_comparison(detected_centers, detected_radii, design_centers, design_radii, matched_pairs, missing_design, extra_detected)
% 可视化比较结果

figure('Name', '检测圆孔与设计圆孔比较');
hold on;

% 计算显示范围
all_centers = [detected_centers; design_centers];
all_radii = [detected_radii; design_radii];
min_x = min(all_centers(:,1) - all_radii);
max_x = max(all_centers(:,1) + all_radii);
min_y = min(all_centers(:,2) - all_radii);
max_y = max(all_centers(:,2) + all_radii);

% 绘制所有设计圆孔 (蓝色实线)
theta = linspace(0, 2*pi, 100);
for i = 1:size(design_centers, 1)
    x = design_centers(i,1) + design_radii(i) * cos(theta);
    y = design_centers(i,2) + design_radii(i) * sin(theta);
    plot(x, y, 'b-', 'LineWidth', 1.5);
end

% 绘制所有检测到的圆孔 (红色虚线)
for i = 1:size(detected_centers, 1)
    x = detected_centers(i,1) + detected_radii(i) * cos(theta);
    y = detected_centers(i,2) + detected_radii(i) * sin(theta);
    plot(x, y, 'r--', 'LineWidth', 1.0);
end

% 高亮显示匹配的圆对
for i = 1:size(matched_pairs, 1)
    det_idx = matched_pairs(i, 1);
    des_idx = matched_pairs(i, 2);
    
    % 连接匹配对的中心点
    plot([detected_centers(det_idx,1), design_centers(des_idx,1)], ...
         [detected_centers(det_idx,2), design_centers(des_idx,2)], ...
         'g-', 'LineWidth', 0.5);
    
    % 在匹配对中心标记编号
    text_pos_x = (detected_centers(det_idx,1) + design_centers(des_idx,1))/2;
    text_pos_y = (detected_centers(det_idx,2) + design_centers(des_idx,2))/2;
    text(text_pos_x, text_pos_y, num2str(i), 'Color', 'k', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

% 标记未被检测到的设计圆孔 (黄色X)
for i = 1:length(missing_design)
    idx = missing_design(i);
    plot(design_centers(idx,1), design_centers(idx,2), 'yx', 'MarkerSize', 10, 'LineWidth', 2);
end

% 标记多余的检测圆孔 (紫色三角)
for i = 1:length(extra_detected)
    idx = extra_detected(i);
    plot(detected_centers(idx,1), detected_centers(idx,2), 'm^', 'MarkerSize', 8, 'LineWidth', 2);
end

% 设置坐标轴和标题
axis equal;
grid on;
xlim([min_x - 10, max_x + 10]);
ylim([min_y - 10, max_y + 10]);
xlabel('X坐标 (设计单位)');
ylabel('Y坐标 (设计单位)');
title('检测圆孔与设计圆孔比较');

% 添加图例
legend('设计圆孔', '检测圆孔', '匹配连线', '未检测到的设计圆孔', '多余检测的圆孔', ...
    'Location', 'best');

hold off;

% 添加误差统计图
if ~isempty(matched_pairs)
    figure('Name', '匹配圆孔的误差统计');
    
    % 中心位置误差直方图
    subplot(2,1,1);
    histogram(sqrt(sum((detected_centers(matched_pairs(:,1),:) - design_centers(matched_pairs(:,2),:)).^2, 2)));
    title('中心位置误差直方图');
    xlabel('中心位置误差 (设计单位)');
    ylabel('频率');
    grid on;
    
    % 半径误差直方图
    subplot(2,1,2);
    histogram(detected_radii(matched_pairs(:,1)) - design_radii(matched_pairs(:,2)));
    title('半径误差直方图');
    xlabel('半径误差 (设计单位)');
    ylabel('频率');
    grid on;
end
end
