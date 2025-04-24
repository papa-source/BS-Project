function [centers, radii] = extract_circle_simple(dxf_file)
% EXTRACT_CIRCLE_SIMPLE 简单的DXF圆提取函数，专为R12格式设计
%
% 输入参数:
%   dxf_file - DXF文件路径，如果未提供，将弹出文件选择对话框
%
% 输出:
%   centers - Nx2矩阵，包含所有圆的圆心坐标 [X,Y]
%   radii - Nx1向量，包含所有圆的半径

% 初始化输出
centers = [];
radii = [];

% 如果未提供文件路径，弹出文件选择对话框
if nargin < 1
    [filename, pathname] = uigetfile({'*.dxf', 'DXF文件 (*.dxf)'}, '选择DXF文件');
    if isequal(filename, 0) || isequal(pathname, 0)
        disp('用户取消了操作');
        return;
    end
    dxf_file = fullfile(pathname, filename);
    disp(['已选择DXF文件: ', dxf_file]);
end

% 读取DXF文件
try
    fid = fopen(dxf_file, 'r');
    if fid == -1
        error('无法打开文件');
    end
    
    % 读取所有内容
    all_lines = textscan(fid, '%s', 'Delimiter', '\n');
    all_lines = all_lines{1};
    fclose(fid);
    
    disp(['文件包含 ', num2str(length(all_lines)), ' 行']);
catch
    disp('读取文件时出错');
    if fid ~= -1
        fclose(fid);
    end
    return;
end

% 极其简单的解析方法
disp('使用超简单的解析方法...');
i = 1;
circle_count = 0;

while i < length(all_lines)
    % 如果找到CIRCLE实体
    if strcmp(strtrim(all_lines{i}), 'CIRCLE')
        disp(['在第 ', num2str(i), ' 行找到CIRCLE']);
        
        % 查找圆参数
        center_x = NaN;
        center_y = NaN;
        radius = NaN;
        
        % 搜索接下来的行
        j = i;
        search_complete = false;
        
        while j < min(i+50, length(all_lines)) && ~search_complete
            j = j + 1;
            
            % 如果找到新实体，停止搜索
            if strcmp(strtrim(all_lines{j}), '0') && j > i+3
                search_complete = true;
                continue;
            end
            
            % 圆心X坐标（组码10）
            if strcmp(strtrim(all_lines{j}), '10')
                if j+1 <= length(all_lines)
                    center_x = str2double(strtrim(all_lines{j+1}));
                    disp(['  找到圆心X: ', num2str(center_x)]);
                end
            end
            
            % 圆心Y坐标（组码20）
            if strcmp(strtrim(all_lines{j}), '20')
                if j+1 <= length(all_lines)
                    center_y = str2double(strtrim(all_lines{j+1}));
                    disp(['  找到圆心Y: ', num2str(center_y)]);
                end
            end
            
            % 半径（组码40）
            if strcmp(strtrim(all_lines{j}), '40')
                if j+1 <= length(all_lines)
                    radius = str2double(strtrim(all_lines{j+1}));
                    disp(['  找到半径: ', num2str(radius)]);
                end
            end
        end
        
        % 如果找到了完整的圆参数，添加到结果
        if ~isnan(center_x) && ~isnan(center_y) && ~isnan(radius)
            centers = [centers; [center_x, center_y]];
            radii = [radii; radius];
            circle_count = circle_count + 1;
            disp(['  成功添加圆 #', num2str(circle_count), ': 圆心=(', num2str(center_x), ', ', num2str(center_y), '), 半径=', num2str(radius)]);
        else
            disp('  未找到完整的圆参数，跳过此圆');
        end
    end
    
    i = i + 1;
end

% 显示结果
disp(['共找到 ', num2str(circle_count), ' 个圆']);

% 如果找到了圆，显示可视化结果
if ~isempty(centers)
    visualize_circles(centers, radii);
end

end

function visualize_circles(centers, radii)
% 可视化找到的圆
figure('Name', 'DXF文件中的圆');
hold on;

% 计算显示范围
min_x = min(centers(:,1) - radii);
max_x = max(centers(:,1) + radii);
min_y = min(centers(:,2) - radii);
max_y = max(centers(:,2) + radii);

% 绘制每个圆
for i = 1:size(centers, 1)
    % 生成圆周上的点
    theta = linspace(0, 2*pi, 100);
    x = centers(i,1) + radii(i) * cos(theta);
    y = centers(i,2) + radii(i) * sin(theta);
    
    % 绘制圆
    plot(x, y, 'b-', 'LineWidth', 1.5);
    
    % 在圆心处标记
    plot(centers(i,1), centers(i,2), 'r+', 'MarkerSize', 10);
    
    % 标记圆的编号
    text(centers(i,1), centers(i,2), num2str(i), 'Color', 'k', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

% 设置坐标轴属性
axis equal;
grid on;

% 考虑边界情况，确保不会出错
if min_x < max_x && min_y < max_y
    margin = max(max_x - min_x, max_y - min_y) * 0.1; % 10%的边距
    xlim([min_x - margin, max_x + margin]);
    ylim([min_y - margin, max_y + margin]);
else
    xlim([-10, 10]);
    ylim([-10, 10]);
end

xlabel('X坐标');
ylabel('Y坐标');
title('DXF文件中的圆');

% 添加图例
legend('圆轮廓', '圆心', 'Location', 'best');

hold off;
end
