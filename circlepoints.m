function [x, y] = circlepoints(cx, cy, r)
% CIRCLEPOINTS 生成圆周上的点坐标
%   [X, Y] = CIRCLEPOINTS(CX, CY, R) 返回中心点为(CX,CY)，半径为R的圆周上的点坐标
%   使用Bresenham算法产生圆周上的点

theta = linspace(0, 2*pi, ceil(2*pi*r));
x = round(cx + r * cos(theta));
y = round(cy + r * sin(theta));
end
