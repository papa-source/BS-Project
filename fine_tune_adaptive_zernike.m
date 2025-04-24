% 细化自适应Zernike算法参数优化
% 在最佳点附近进行更精细的参数搜索

% 加载边缘点数据
load artificial_moment_data.mat ex1 ey1 ex2 ey2 theory_c1 theory_c2

% 使用理论值作为参考
true_c1 = theory_c1;
true_c2 = theory_c2;

% 初始化更精细的参数范围 - 围绕最佳点（相位=0, 幅度因子=2.0）
phase_offsets = [-0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2];  % 微小相位偏移
magnitude_factors = [1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2];  % 幅度因子
radius_factors = [0, 0.01, 0.02, 0.03, 0.04, 0.05];  % 半径调整因子

% 测试更多不同的归一化方法
normalization_methods = {'standard', 'robust', 'adaptive'};
current_method = 'standard';

% 存储结果
num_configs = length(phase_offsets) * length(magnitude_factors) * length(radius_factors);
results = zeros(num_configs, 5);  % [phase_idx, mag_idx, rad_idx, 误差1, 误差2]
result_idx = 1;

fprintf('细化参数优化开始，共%d种配置...\n', num_configs);

% 遍历所有参数组合
for p_idx = 1:length(phase_offsets)
    phase_offset = phase_offsets(p_idx);
    
    for m_idx = 1:length(magnitude_factors)
        mag_factor = magnitude_factors(m_idx);
        
        for r_idx = 1:length(radius_factors)
            rad_factor = radius_factors(r_idx);
            
            % 测试当前参数组合
            [c1, c2] = test_config(ex1, ey1, ex2, ey2, phase_offset, mag_factor, rad_factor, current_method);
            
            % 计算误差
            error1 = sqrt(sum((c1 - true_c1).^2));
            error2 = sqrt(sum((c2 - true_c2).^2));
            avg_error = (error1 + error2) / 2;
            
            % 存储结果
            results(result_idx, :) = [p_idx, m_idx, r_idx, error1, error2];
            
            fprintf('配置 %d/%d: 相位偏移=%.2f°, 幅度因子=%.2f, 半径因子=%.3f, 平均误差=%.4f像素\n', ...
                result_idx, num_configs, phase_offset*180/pi, mag_factor, rad_factor, avg_error);
            
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
fprintf('半径因子 = %.3f\n', best_rad);
fprintf('平均误差 = %.4f像素\n', mean(best_config(4:5)));
fprintf('误差1 = %.4f像素, 误差2 = %.4f像素\n', best_config(4), best_config(5));

% 保存最佳配置
save('best_fine_tuned_params.mat', 'best_phase', 'best_mag', 'best_rad');

% 可视化结果
visualize_results(results, phase_offsets, magnitude_factors, radius_factors);

% 测试不同归一化方法
fprintf('\n\n测试不同归一化方法的影响...\n');
for i = 1:length(normalization_methods)
    method = normalization_methods{i};
    [c1, c2] = test_config(ex1, ey1, ex2, ey2, best_phase, best_mag, best_rad, method);
    error1 = sqrt(sum((c1 - true_c1).^2));
    error2 = sqrt(sum((c2 - true_c2).^2));
    avg_error = (error1 + error2) / 2;
    fprintf('归一化方法: %s, 平均误差: %.4f像素\n', method, avg_error);
end

%% 辅助函数

% 使用参数配置测试自适应Zernike
function [c1, c2] = test_config(ex1, ey1, ex2, ey2, phase_offset, mag_factor, rad_factor, norm_method)
    % 第一个圆
    center_est1 = [mean(ex1), mean(ey1)];
    radius_est1 = mean(sqrt((ex1-center_est1(1)).^2 + (ey1-center_est1(2)).^2));
    
    % 根据归一化方法选择不同的处理方式
    if strcmp(norm_method, 'standard')
        % 标准归一化
        norm_points1 = [(ex1 - center_est1(1))/radius_est1, (ey1 - center_est1(2))/radius_est1];
    elseif strcmp(norm_method, 'robust')
        % 鲁棒归一化 - 使用中位数和四分位距
        median_x = median(ex1);
        median_y = median(ey1);
        iqr_x = iqr(ex1);
        iqr_y = iqr(ey1);
        if iqr_x == 0 || iqr_y == 0
            norm_points1 = [(ex1 - median_x)/radius_est1, (ey1 - median_y)/radius_est1];
        else
            norm_points1 = [(ex1 - median_x)/iqr_x, (ey1 - median_y)/iqr_y];
        end
    else
        % 自适应归一化 - 使用边缘点的方差
        std_x = std(ex1);
        std_y = std(ey1);
        if std_x == 0 || std_y == 0
            norm_points1 = [(ex1 - center_est1(1))/radius_est1, (ey1 - center_est1(2))/radius_est1];
        else
            norm_points1 = [(ex1 - center_est1(1))/std_x, (ey1 - center_est1(2))/std_y];
        end
    end
    
    % 计算Zernike矩并使用改进的实现
    Z31_1 = calculate_zernike_moment(norm_points1, 3, 1);
    
    % 计算圆心修正
    phase_Z31 = angle(Z31_1);
    magnitude_Z31 = abs(Z31_1);
    
    % 应用参数
    dx1 = radius_est1 * magnitude_Z31 * cos(phase_Z31 + phase_offset) / mag_factor;
    dy1 = radius_est1 * magnitude_Z31 * sin(phase_Z31 + phase_offset) / mag_factor;
    
    % 更新圆心坐标和半径
    c1 = center_est1 + [dx1, dy1];
    
    % 第二个圆 - 过程类似
    center_est2 = [mean(ex2), mean(ey2)];
    radius_est2 = mean(sqrt((ex2-center_est2(1)).^2 + (ey2-center_est2(2)).^2));
    
    if strcmp(norm_method, 'standard')
        norm_points2 = [(ex2 - center_est2(1))/radius_est2, (ey2 - center_est2(2))/radius_est2];
    elseif strcmp(norm_method, 'robust')
        median_x = median(ex2);
        median_y = median(ey2);
        iqr_x = iqr(ex2);
        iqr_y = iqr(ey2);
        if iqr_x == 0 || iqr_y == 0
            norm_points2 = [(ex2 - median_x)/radius_est2, (ey2 - median_y)/radius_est2];
        else
            norm_points2 = [(ex2 - median_x)/iqr_x, (ey2 - median_y)/iqr_y];
        end
    else
        std_x = std(ex2);
        std_y = std(ey2);
        if std_x == 0 || std_y == 0
            norm_points2 = [(ex2 - center_est2(1))/radius_est2, (ey2 - center_est2(2))/radius_est2];
        else
            norm_points2 = [(ex2 - center_est2(1))/std_x, (ey2 - center_est2(2))/std_y];
        end
    end
    
    Z31_2 = calculate_zernike_moment(norm_points2, 3, 1);
    
    phase_Z31 = angle(Z31_2);
    magnitude_Z31 = abs(Z31_2);
    
    dx2 = radius_est2 * magnitude_Z31 * cos(phase_Z31 + phase_offset) / mag_factor;
    dy2 = radius_est2 * magnitude_Z31 * sin(phase_Z31 + phase_offset) / mag_factor;
    
    c2 = center_est2 + [dx2, dy2];
end

% 改进的Zernike矩计算
function Z = calculate_zernike_moment(points, n, m)
    % 提取点坐标
    x = points(:, 1);
    y = points(:, 2);
    
    % 转换为极坐标
    rho = sqrt(x.^2 + y.^2);
    theta = atan2(y, x);
    
    % 优化：确保全部使用单位圆内的有效点
    valid_idx = rho <= 1;
    if sum(valid_idx) < 3
        % 如果有效点太少，调整缩放以包含更多点
        max_rho = max(rho);
        if max_rho > 0
            rho = rho / max_rho;
            valid_idx = ones(size(rho));
        else
            Z = 0;
            return;
        end
    else
        % 只使用单位圆内的点
        rho = rho(valid_idx);
        theta = theta(valid_idx);
    end
    
    % 边缘点默认亮度为1
    f = ones(size(rho));
    
    % 计算径向多项式 - 针对n=3,m=1进行特殊优化
    if n == 3 && abs(m) == 1
        % R31 = (3ρ³ - 2ρ) - 直接使用闭式解以提高精度
        R = 3*rho.^3 - 2*rho;
    else
        % 标准计算其他阶数
        R = zeros(size(rho));
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
        Vreal = R .* cos(m*theta);
        Vimag = R .* sin(m*theta);
    else
        Vreal = R .* cos(abs(m)*theta);
        Vimag = -R .* sin(abs(m)*theta);
    end
    
    % 计算Zernike矩 - 使用点数加权以提高稳定性
    Z_real = sum(f .* Vreal) / length(f);
    Z_imag = sum(f .* Vimag) / length(f);
    
    % 组合成复数
    Z = complex(Z_real, Z_imag);
    
    % 确保结果有效
    if isnan(Z) || isinf(Z)
        Z = 0;
    end
end

% 可视化结果
function visualize_results(results, phase_offsets, magnitude_factors, radius_factors)
    figure('Name', '自适应Zernike参数优化', 'Position', [100, 100, 1200, 600]);
    
    % 为每个半径因子创建一个子图
    subplot_rows = 2;
    subplot_cols = ceil(length(radius_factors)/2);
    
    for r_idx = 1:length(radius_factors)
        subplot(subplot_rows, subplot_cols, r_idx);
        
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
        title(sprintf('半径因子 = %.3f', radius_factors(r_idx)));
        
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
    saveas(gcf, 'fine_tuned_zernike_optimization.png');
end
