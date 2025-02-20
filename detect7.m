function [centers, radii] = detect7(edges)
% 圆孔识别与拟合
% 使用边缘检测结果进行圆孔检测和拟合

% 使用边缘检测结果进行圆检测
[centers, radii] = imfindcircles(edges, [10 30], ...  % 调整半径范围为10-30像素
    'ObjectPolarity', 'bright', ...
    'Sensitivity', 0.85, ...  % 降低敏感度
    'EdgeThreshold', 0.1, ...
    'Method', 'PhaseCode');

% 如果没有检测到圆，尝试使用二值掩码
if isempty(centers)
    % 注意：这里需要调用方提供holes_mask参数，或者在函数声明中添加该参数
    % 为了保持兼容性，这里我们返回空结果
    centers = [];
    radii = [];
end
end