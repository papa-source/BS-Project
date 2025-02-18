function [lines, corners] = detect_and_fit_contours(edges)
% 轮廓线识别与拟合
% 输入参数:
%   edges - 边缘图像
% 输出参数:
%   lines - 检测到的直线
%   corners - 检测到的角点

% 1. 使用Hough变换检测直线
[H, theta, rho] = hough(edges);
P = houghpeaks(H, 10, 'threshold', ceil(0.3*max(H(:))));
lines = houghlines(edges, theta, rho, P, 'FillGap', 5, 'MinLength', 7);

% 2. 获取角点
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
corner = [x, y];
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