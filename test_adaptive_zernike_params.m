% 测试自适应Zernike算法的不同参数配置
% 此脚本测试不同的相位偏移、幅度因子和半径因子，找出最佳组合

% 加载已知的真实圆心位置（从人工生成的测试数据）
true_c1 = [190.11, 55.89];
true_c2 = [90.45, 155.34];

% 加载预处理好的边缘点数据
load artificial_moment_data.mat ex1 ey1 ex2 ey2 theory_c1 theory_c2

% 使用真实的圆心作为参考
true_c1 = theory_c1;
true_c2 = theory_c2;

% 初始化参数范围
phase_offsets = [0, pi/6, pi/4, pi/3, pi/2, 2*pi/3, 3*pi/4, 5*pi/6, pi];  % 相位偏移
magnitude_factors = [2, 2.5, 3, 3.5, 4];  % 幅度因子
radius_factors = [0, 0.1, 0.15, 0.2, 0.25, 0.3];  % 半径调整因子

% 存储结果
num_configs = length(phase_offsets) * length(magnitude_factors) * length(radius_factors);
results = zeros(num_configs, 5);  % [phase_idx, mag_idx, rad_idx, 误差1, 误差2]
result_idx = 1;

fprintf('参数优化测试开始，共%d种配置...\n', num_configs);

% 遍历所有参数组合
for p_idx = 1:length(phase_offsets)
    phase_offset = phase_offsets(p_idx);
    
    for m_idx = 1:length(magnitude_factors)
        mag_factor = magnitude_factors(m_idx);
        
        for r_idx = 1:length(radius_factors)
            rad_factor = radius_factors(r_idx);
            
            % 测试当前参数组合
            [c1, c2] = test_single_config(ex1, ey1, ex2, ey2, phase_offset, mag_factor, rad_factor);
            
            % 计算误差
            error1 = sqrt(sum((c1 - true_c1).^2));
            error2 = sqrt(sum((c2 - true_c2).^2));
            avg_error = (error1 + error2) / 2;
            
            % 存储结果
            results(result_idx, :) = [p_idx, m_idx, r_idx, error1, error2];
            
            fprintf('配置 %d/%d: 相位偏移=%.2f (%.1f°), 幅度因子=%.1f, 半径因子=%.2f, 平均误差=%.4f像素\n', ...
                result_idx, num_configs, phase_offset, phase_offset*180/pi, mag_factor, rad_factor, avg_error);
            
            result_idx = result_idx + 1;
        end
    end
end

% 找出最佳配置
[~, best_idx] = min(mean(results(:, 4:5), 2));
best_config = results(best_idx, :);
best_phase = phase_offsets(best_config(1));
best_mag = magnitude_factors(best_config(2));
best_rad = radius_factors(best_config(3));

fprintf('\n最佳参数配置:\n');
fprintf('相位偏移 = %.4f (%.2f°)\n', best_phase, best_phase*180/pi);
fprintf('幅度因子 = %.2f\n', best_mag);
fprintf('半径因子 = %.2f\n', best_rad);
fprintf('平均误差 = %.4f像素\n', mean(best_config(4:5)));
fprintf('误差1 = %.4f像素, 误差2 = %.4f像素\n', best_config(4), best_config(5));

% 保存最佳配置
save('best_adaptive_zernike_params.mat', 'best_phase', 'best_mag', 'best_rad');

% 绘制误差热图
visualize_results(results, phase_offsets, magnitude_factors, radius_factors);

%% 辅助函数

% 使用单一参数配置测试自适应Zernike
function [c1, c2] = test_single_config(ex1, ey1, ex2, ey2, phase_offset, mag_factor, rad_factor)
    % 第一个圆
    center_est1 = [mean(ex1), mean(ey1)];
    radius_est1 = mean(sqrt((ex1-center_est1(1)).^2 + (ey1-center_est1(2)).^2));
    norm_points1 = [(ex1 - center_est1(1))/radius_est1, (ey1 - center_est1(2))/radius_est1];
    Z31_1 = calculate_zernike_moment(norm_points1, 3, 1);
    
    % 使用传入的参数计算圆心修正
    phase_Z31 = angle(Z31_1);
    magnitude_Z31 = abs(Z31_1);
    
    % 应用参数
    dx1 = radius_est1 * magnitude_Z31 * cos(phase_Z31 + phase_offset) / mag_factor;
    dy1 = radius_est1 * magnitude_Z31 * sin(phase_Z31 + phase_offset) / mag_factor;
    
    % 更新圆心坐标
    c1 = center_est1 + [dx1, dy1];
    
    % 第二个圆
    center_est2 = [mean(ex2), mean(ey2)];
    radius_est2 = mean(sqrt((ex2-center_est2(1)).^2 + (ey2-center_est2(2)).^2));
    norm_points2 = [(ex2 - center_est2(1))/radius_est2, (ey2 - center_est2(2))/radius_est2];
    Z31_2 = calculate_zernike_moment(norm_points2, 3, 1);
    
    % 使用传入的参数计算圆心修正
    phase_Z31 = angle(Z31_2);
    magnitude_Z31 = abs(Z31_2);
    
    % 应用参数
    dx2 = radius_est2 * magnitude_Z31 * cos(phase_Z31 + phase_offset) / mag_factor;
    dy2 = radius_est2 * magnitude_Z31 * sin(phase_Z31 + phase_offset) / mag_factor;
    
    % 更新圆心坐标
    c2 = center_est2 + [dx2, dy2];
end

% 可视化参数测试结果
function visualize_results(results, phase_offsets, magnitude_factors, radius_factors)
    figure('Name', '自适应Zernike参数优化', 'Position', [100, 100, 1200, 600]);
    
    % 为每个半径因子创建一个子图
    num_rad = length(radius_factors);
    for r_idx = 1:num_rad
        subplot(2, ceil(num_rad/2), r_idx);
        
        % 提取当前半径因子的数据
        mask = results(:, 3) == r_idx;
        current_results = results(mask, :);
        
        % 创建热图数据
        error_map = zeros(length(phase_offsets), length(magnitude_factors));
        for i = 1:size(current_results, 1)
            p_idx = current_results(i, 1);
            m_idx = current_results(i, 2);
            error_map(p_idx, m_idx) = mean(current_results(i, 4:5));
        end
        
        % 绘制热图
        imagesc(error_map);
        colormap(jet);
        colorbar;
        
        % 设置标签
        xlabel('幅度因子');
        ylabel('相位偏移');
        title(sprintf('半径因子 = %.2f', radius_factors(r_idx)));
        
        % 设置刻度标签
        xticks(1:length(magnitude_factors));
        xticklabels(arrayfun(@num2str, magnitude_factors, 'UniformOutput', false));
        yticks(1:length(phase_offsets));
        yticklabels(arrayfun(@(x) sprintf('%.1f°', x*180/pi), phase_offsets, 'UniformOutput', false));
        
        % 标记最小误差点
        [min_val, min_idx] = min(reshape(error_map, [], 1));
        [min_y, min_x] = ind2sub(size(error_map), min_idx);
        hold on;
        plot(min_x, min_y, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        text(min_x, min_y, sprintf('  %.4f', min_val), 'Color', 'r', 'FontWeight', 'bold');
        hold off;
    end
    
    sgtitle('自适应Zernike参数优化 - 平均误差 (像素)');
    saveas(gcf, 'adaptive_zernike_parameter_optimization.png');
end

% 计算Zernike矩的辅助函数
function Z = calculate_zernike_moment(points, n, m)
    % 计算Zernike矩 - 实现论文中的公式(4-23)(4-32)
    % Z_nm = ∑∑f(x,y)·V*_nm(ρ,θ)
    
    % 提取点坐标
    x = points(:, 1);
    y = points(:, 2);
    
    % 转换为极坐标
    rho = sqrt(x.^2 + y.^2);
    theta = atan2(y, x);
    
    % 仅使用单位圆内的点
    valid_idx = rho <= 1;
    if sum(valid_idx) < 5  % 至少需要5个有效点
        Z = 0;
        return;
    end
    
    rho = rho(valid_idx);
    theta = theta(valid_idx);
    
    % 边缘点默认亮度为1
    f = ones(size(rho));
    
    % 计算径向多项式
    R = zeros(size(rho));
    
    % 根据阶数n和角度m计算径向多项式
    % 特殊情况处理
    if n == 0 && m == 0
        % R00 = 1
        R = ones(size(rho));
    elseif n == 1 && abs(m) == 1
        % R11 = ρ
        R = rho;
    elseif n == 2 && m == 0
        % R20 = 2ρ² - 1
        R = 2*rho.^2 - 1;
    elseif n == 2 && abs(m) == 2
        % R22 = ρ²
        R = rho.^2;
    elseif n == 3 && abs(m) == 1
        % R31 = (3ρ³ - 2ρ)
        R = 3*rho.^3 - 2*rho;
    elseif n == 3 && abs(m) == 3
        % R33 = ρ³
        R = rho.^3;
    else
        % 一般情况的径向多项式
        for s = 0:((n-abs(m))/2)
            coef = (-1)^s * factorial(n-s) / ...
                  (factorial(s) * factorial((n+abs(m))/2-s) * factorial((n-abs(m))/2-s));
            
            % 防止数值溢出
            if isnan(coef) || isinf(coef)
                continue;
            end
            
            R = R + coef * rho.^(n-2*s);
        end
    end
    
    % 计算复值Zernike矩
    if m >= 0
        % Vnm = Rnm(ρ)·exp(jmθ) = Rnm(ρ)·(cos(mθ) + j·sin(mθ))
        Vreal = R .* cos(m*theta);
        Vimag = R .* sin(m*theta);
    else
        % Vn,-m = Rnm(ρ)·exp(-jmθ) = Rnm(ρ)·(cos(mθ) - j·sin(mθ))
        Vreal = R .* cos(abs(m)*theta);
        Vimag = -R .* sin(abs(m)*theta);
    end
    
    % 计算Zernike矩
    Z_real = sum(f .* Vreal) / sum(valid_idx);
    Z_imag = sum(f .* Vimag) / sum(valid_idx);
    
    % 组合成复数
    Z = complex(Z_real, Z_imag);
    
    % 确保结果有效
    if isnan(Z) || isinf(Z)
        Z = 0;
        warning('Zernike矩计算结果是NaN或Inf，使用默认值0');
    end
end
