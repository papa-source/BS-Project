function display_results(original_img, edges, centers, radii, lines, corners)
% 显示检测结果
figure('Name', '板材冲孔质量检测结果', 'NumberTitle', 'off');

% 1. 原图
subplot(2,2,1);
imshow(original_img);
title('原始图像', 'FontSize', 12);

% 2. 边缘检测结果
subplot(2,2,2);
imshow(edges);
title('边缘检测结果', 'FontSize', 12);

% 3. 圆孔检测结果
subplot(2,2,3);
imshow(original_img);
hold on;
if ~isempty(centers)
    % 先画圆
    viscircles(centers, radii, 'EdgeColor', 'b', 'LineWidth', 1.5);
    
    % 对圆孔进行排序（从左上到右下）
    [~, order] = sortrows(centers, [2 1]);  % 先按y排序，再按x排序
    
    % 标注圆心和半径，使用更清晰的布局
    for i = 1:size(centers, 1)
        idx = order(i);
        % 画圆心
        plot(centers(idx,1), centers(idx,2), 'b+', 'MarkerSize', 8, 'LineWidth', 1.5);
        
        % 计算标注位置（避免重叠）
        angle = 45;  % 标注位置的角度
        offset = radii(idx) * 1.2;  % 标注距离圆心的距离
        text_x = centers(idx,1) + offset * cosd(angle);
        text_y = centers(idx,2) - offset * sind(angle);
        
        % 添加标注文本
        text(text_x, text_y, sprintf('%d\nR=%.1f', i, radii(idx)), ...
             'Color', 'blue', 'FontSize', 8, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'left', ...
             'BackgroundColor', [1 1 1 0.7]);  % 半透明白色背景
    end
end
title('圆孔检测结果', 'FontSize', 12);
hold off;

% 4. 轮廓线检测结果
subplot(2,2,4);
imshow(original_img);
hold on;

% 绘制轮廓线
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1), xy(:,2), 'LineWidth', 2, 'Color', [0 0.7 0]);
end

% 绘制角点
if ~isempty(corners)
    plot(corners(:,1), corners(:,2), 'r*', 'MarkerSize', 8);
    % 标注角点编号
    for i = 1:size(corners, 1)
        text(corners(i,1)+5, corners(i,2)+5, sprintf('%d', i), ...
             'Color', 'red', 'FontSize', 8, 'FontWeight', 'bold', ...
             'BackgroundColor', [1 1 1 0.7]);
    end
end

title('轮廓线检测结果', 'FontSize', 12);
hold off;

% 调整图像显示
set(gcf, 'Position', get(0, 'Screensize'));  % 全屏显示
end