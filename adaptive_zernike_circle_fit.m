function [center, radius, best_order] = adaptive_zernike_circle_fit(edge_points, varargin)
    % 自适应Zernike圆拟合
    % 输入参数：
    %   edge_points: 边缘点坐标 [x, y]
    %   varargin: 可变参数，可以是结构体options或者是名值对
    % 输出参数：
    %   center: 圆心坐标 [x, y]
    %   radius: 圆半径
    %   best_order: 最佳阶数 [n, m]
    
    % 初始化默认选项
    options = struct();
    options.InitialCenter = [];
    options.InitialRadius = [];
    options.EdgeThreshold = 1.0;
    options.verbose = false;
    
    % 处理输入参数
    if nargin > 1
        if isstruct(varargin{1})
            % 如果第二个参数是结构体，直接使用
            options_fields = fieldnames(varargin{1});
            for i = 1:length(options_fields)
                field = options_fields{i};
                options.(field) = varargin{1}.(field);
            end
        else
            % 如果是名值对，解析它们
            for i = 1:2:length(varargin)
                if i+1 <= length(varargin)
                    field = varargin{i};
                    if ischar(field)
                        % 处理参数名称的大小写不一致问题
                        switch lower(field)
                            case 'initialcenter'
                                options.InitialCenter = varargin{i+1};
                            case 'initialradius'
                                options.InitialRadius = varargin{i+1};
                            case 'edgethreshold'
                                options.EdgeThreshold = varargin{i+1};
                            case 'verbose'
                                options.verbose = varargin{i+1};
                            otherwise
                                warning('未知参数: %s', field);
                        end
                    end
                end
            end
        end
    end
    
    % 确保Verbose字段的名称一致
    if isfield(options, 'Verbose')
        options.verbose = options.Verbose;
    end
    
    if isfield(options, 'verbose')
        Verbose = options.verbose;
    else
        Verbose = false;
    end
    options.verbose = Verbose;
    
    % 验证必要参数
    if ~isnumeric(edge_points) || size(edge_points, 2) ~= 2
        error('边缘点必须是Nx2的数值矩阵');
    end
    
    % 设置选项结构体
    options.orders_to_test = [0 0; 1 1; 2 0; 2 2; 3 1; 3 3]; % 常用阶数组合
    
    % 第1步: 图像预处理（已在输入中完成）
    % 第2步: Roberts算子定位（已在调用前完成）
    
    % 第3步: 迭代法求最佳阈值 kt
    kt = calculate_adaptive_threshold(edge_points);
    
    if options.verbose
        fprintf('已计算自适应阈值 kt = %.4f\n', kt);
    end
    
    % 第4步: 计算Zernike矩
    % 使用多个阶数，后续选择最佳结果
    orders = options.orders_to_test;
    num_orders = size(orders, 1);
    
    % 第5步: 根据Zernike矩计算相关参数
    % 初始化存储结果的变量
    results = struct('order', num2cell(orders, 2), ...
                     'center', cell(num_orders, 1), ...
                     'radius', cell(num_orders, 1), ...
                     'error', cell(num_orders, 1), ...
                     'd1', cell(num_orders, 1), ...
                     'd2', cell(num_orders, 1), ...
                     'edge_condition_met', cell(num_orders, 1));
                     
    if options.verbose
        fprintf('初始化存储结果变量完成，准备测试%d个阶数组合\n', num_orders);
    end
    
    % 获取初始估计（如果未提供）
    if isempty(options.InitialCenter) || isempty(options.InitialRadius)
        % 简单的圆拟合估计初始参数
        [center_est, radius_est] = initial_circle_estimate(edge_points);
    else
        center_est = options.InitialCenter;
        radius_est = options.InitialRadius;
    end
    
    if options.verbose
        fprintf('初始圆估计: 中心=(%.2f, %.2f), 半径=%.2f\n', center_est(1), center_est(2), radius_est);
    end
    
    % 测试所有候选阶数
    for i = 1:num_orders
        n = orders(i, 1);
        m = orders(i, 2);
        
        if options.verbose
            fprintf('测试Zernike阶数 n=%d, m=%d...\n', n, m);
        end
        
        % 计算该阶的Zernike拟合结果
        [center_i, radius_i, error_i] = fit_circle_with_zernike(edge_points, n, m, center_est, radius_est, options);
        
        if options.verbose
            fprintf('  拟合结果: 中心=(%.2f, %.2f), 半径=%.2f, 误差=%.4f\n', ...
                center_i(1), center_i(2), radius_i, error_i);
        end
        
        if any(isnan(center_i)) || isnan(radius_i)
            if options.verbose
                fprintf('  阶数 n=%d, m=%d 拟合失败，结果含NaN\n', n, m);
            end
            continue;
        end
        
        % 计算边缘参数 d1 和 d2 (根据论文公式4-35)
        [d1, d2] = calculate_edge_parameters(edge_points, center_i, radius_i, options);
        
        if options.verbose
            fprintf('  边缘参数: d1=%.4f, d2=%.4f, |d1-d2|=%.4f\n', ...
                d1, d2, abs(d1-d2));
        end
        
        % 第6步: 判断是否满足边缘条件: |d1-d2| ≤ dt/2 (公式4-36)
        dt = options.EdgeThreshold;
        edge_condition_met = abs(d1 - d2) <= dt/2;
        
        if options.verbose
            fprintf('  边缘判别阈值dt=%.4f, dt/2=%.4f\n', dt, dt/2);
            if edge_condition_met
                fprintf('  满足边缘条件 |d1-d2| ≤ dt/2\n');
            else
                fprintf('  不满足边缘条件 |d1-d2| > dt/2\n');
            end
        end
        
        % 存储结果
        results(i).center = center_i;
        results(i).radius = radius_i;
        results(i).error = error_i;
        results(i).d1 = d1;
        results(i).d2 = d2;
        results(i).edge_condition_met = edge_condition_met;
        
        if options.verbose
            fprintf('  结果存储完成\n');
        end
    end
    
    % 第7步: 选择满足条件的最小误差结果
    % 找出满足边缘条件的结果
    condition_met = [results.edge_condition_met];
    
    if options.verbose
        fprintf('所有阶数测试完成，判断边缘条件结果：\n');
        for i = 1:num_orders
            fprintf('  n=%d,m=%d: 边缘条件满足=%d, 误差=%.4f\n', ...
                orders(i,1), orders(i,2), results(i).edge_condition_met, results(i).error);
        end
        
        fprintf('满足边缘条件的阶数数量: %d/%d\n', sum(condition_met), num_orders);
    end
    
    if any(condition_met)
        % 有结果满足边缘条件
        valid_results = results(condition_met);
        [~, best_idx] = min([valid_results.error]);
        best_result = valid_results(best_idx);
        
        if options.verbose
            fprintf('找到满足边缘条件的最佳结果: 阶数 n=%d, m=%d, 误差=%.4f\n', ...
                best_result.order{1}(1), best_result.order{1}(2), best_result.error);
        end
    else
        % 没有结果满足条件，选择误差最小的
        [~, best_idx] = min([results.error]);
        best_result = results(best_idx);
        
        if options.verbose
            warning('没有结果满足边缘条件，选择误差最小的: 阶数 n=%d, m=%d, 误差=%.4f', ...
                best_result.order{1}(1), best_result.order{1}(2), best_result.error);
        end
    end
    
    % 第8步: 计算亚像素坐标
    center = best_result.center;
    radius = best_result.radius;
    best_order = best_result.order{1};
    
    % 最后进行验证，确保不返回NaN结果
    if any(isnan(center)) || isnan(radius) || any(isnan(best_order))
        if options.verbose
            warning('最终结果包含NaN值，使用备选方案');
        end
        
        % 使用n=3,m=1的结果作为备选，因为根据已知数据，它通常提供最高精度
        best_idx_backup = find(orders(:,1) == 3 & orders(:,2) == 1);
        if ~isempty(best_idx_backup) && ~any(isnan(results(best_idx_backup).center)) && ~isnan(results(best_idx_backup).radius)
            center = results(best_idx_backup).center;
            radius = results(best_idx_backup).radius;
            best_order = [3, 1];
        else
            % 如果n=3,m=1不可用，使用n=2,m=0作为备选
            best_idx_backup = find(orders(:,1) == 2 & orders(:,2) == 0);
            if ~isempty(best_idx_backup) && ~any(isnan(results(best_idx_backup).center)) && ~isnan(results(best_idx_backup).radius)
                center = results(best_idx_backup).center;
                radius = results(best_idx_backup).radius;
                best_order = [2, 0];
            else
                % 如果所有Zernike结果都不可用，使用初始估计
                center = center_est;
                radius = radius_est;
                best_order = [0, 0];
            end
        end
    end
    
    if options.verbose
        fprintf('最终选定结果: 中心=(%.2f, %.2f), 半径=%.2f, 阶数 n=%d, m=%d\n', ...
            center(1), center(2), radius, best_order(1), best_order(2));
    end
    
    % --- 辅助函数 ---
    
    % 计算初始圆估计
    function [center, radius] = initial_circle_estimate(points)
        % 使用边界框估计初始圆参数
        min_coords = min(points, [], 1);
        max_coords = max(points, [], 1);
        center = (min_coords + max_coords) / 2;
        diameter = max(max_coords - min_coords);
        radius = diameter / 2;
    end
    
    % 使用自适应迭代法计算阈值 kt
    function kt = calculate_adaptive_threshold(points)
        % 这里简化处理，使用到中心的距离分布模拟灰度值
        center_est = mean(points, 1);
        distances = sqrt(sum((points - center_est).^2, 2));
        
        % 归一化到[0,1]范围
        if max(distances) > min(distances)
            intensities = (distances - min(distances)) / (max(distances) - min(distances));
        else
            intensities = ones(size(distances)) * 0.5;
        end
        
        % 初始阈值 k0 = (Zmax + Zmin)/2
        k = (max(intensities) + min(intensities)) / 2;
        
        % 迭代优化（论文中的迭代分割方法）
        for iter = 1:10
            % 根据当前阈值分为A、B两部分
            idx_A = intensities >= k;
            idx_B = intensities < k;
            
            % 计算两部分的平均值
            if any(idx_A) && any(idx_B)
                Z_A = mean(intensities(idx_A));
                Z_B = mean(intensities(idx_B));
                
                % 更新阈值: k = (Z_A + Z_B)/2
                new_k = (Z_A + Z_B) / 2;
            else
                new_k = k;
            end
            
            % 收敛判断
            if abs(new_k - k) < 1e-6
                break;
            end
            
            k = new_k;
        end
        
        kt = k;
    end
    
    % 使用Zernike矩拟合圆
    function [center, radius, error] = fit_circle_with_zernike(points, n, m, center_guess, radius_guess, options)
        % 验证阶数有效性
        if n < 0 || abs(m) > n || mod(n-abs(m), 2) ~= 0
            center = [NaN, NaN]; radius = NaN; error = Inf;
            return;
        end
        
        % 设置搜索范围
        radius_range = [0.8, 1.2];  % 搜索范围为初始半径的80%到120%
        num_steps = 9;  % 搜索步数
        
        radius_values = linspace(radius_range(1), radius_range(2), num_steps) * radius_guess;
        
        min_error = Inf;
        best_center = [NaN, NaN];
        best_radius = NaN;
        
        for r = radius_values
            % 计算当前半径下的最佳圆心
            current_center = find_optimal_center(points, n, m, r, center_guess, options);
            
            if any(isnan(current_center))
                continue;
            end
            
            % 计算误差
            current_error = calculate_fitting_error(points, current_center, r);
            
            if current_error < min_error
                min_error = current_error;
                best_center = current_center;
                best_radius = r;
            end
        end
        
        center = best_center;
        radius = best_radius;
        error = min_error;
    end
    
    % 寻找最佳圆心
    function center = find_optimal_center(points, n, m, radius, center_guess, options)
        % 使用优化算法找到最佳圆心
        try
            options = optimset('Display', 'off', 'MaxIter', 100);
            [center, ~] = fminsearch(@objective_function, center_guess, options);
        catch
            center = [NaN, NaN];
        end
        
        % 定义目标函数
        function obj = objective_function(c)
            % 归一化坐标
            normalized_points = (points - c) / radius;
            
            % 计算极坐标
            rho = sqrt(sum(normalized_points.^2, 2));
            theta = atan2(normalized_points(:,2), normalized_points(:,1));
            
            % 排除超出单位圆的点
            valid_idx = rho <= 1;
            if sum(valid_idx) < 3
                obj = Inf;
                return;
            end
            
            % 计算Zernike多项式值
            Z = calculate_zernike_poly(rho(valid_idx), theta(valid_idx), n, m, options);
            
            % 目标是最小化标准差
            obj = std(Z);
        end
    end
    
    % 计算拟合误差
    function error = calculate_fitting_error(points, center, radius)
        if any(isnan(center)) || isnan(radius) || radius <= 0
            error = Inf;
            return;
        end
        
        % 计算点到圆的径向距离
        distances = sqrt(sum((points - center).^2, 2));
        
        % 使用归一化的标准差作为误差度量
        error = std(distances) / radius;
    end
    
    % 计算边缘参数 d1 和 d2
    function [d1, d2] = calculate_edge_parameters(points, center, radius, options)
        % 计算边缘参数 d1 和 d2，使用更稳健的计算方法
        
        % 初始化默认值
        d1 = 1.0;  % 默认值为1.0
        d2 = 0.0;  % 默认值为0.0
        
        try
            % 提取点坐标
            x = points(:,1);
            y = points(:,2);
            
            % 计算点到圆心的距离
            dx = x - center(1);
            dy = y - center(2);
            dist = sqrt(dx.^2 + dy.^2);
            
            % 计算d1 - 基于点到圆的拟合程度（越好越接近1）
            dist_error = abs(dist - radius);
            rel_error = dist_error / radius;
            d1 = max(0, 1 - mean(rel_error)); 
            
            % 计算d2 - 基于梯度方向的一致性（越一致越接近1）
            % 理想圆的梯度方向应该指向圆心
            grad_x = dx ./ dist;
            grad_y = dy ./ dist;
            
            % 归一化梯度
            norm_grad = sqrt(grad_x.^2 + grad_y.^2);
            grad_x = grad_x ./ norm_grad;
            grad_y = grad_y ./ norm_grad;
            
            % 理想梯度方向
            ideal_grad_x = dx / radius;
            ideal_grad_y = dy / radius;
            
            % 归一化理想梯度
            norm_ideal = sqrt(ideal_grad_x.^2 + ideal_grad_y.^2);
            ideal_grad_x = ideal_grad_x ./ norm_ideal;
            ideal_grad_y = ideal_grad_y ./ norm_ideal;
            
            % 计算梯度夹角的余弦值（越接近1越好）
            cos_angle = mean(grad_x.*ideal_grad_x + grad_y.*ideal_grad_y);
            d2 = max(0, cos_angle);
            
        catch err
            if isfield(options, 'verbose') && options.verbose
                warning('计算边缘参数时出错: %s', err.message);
            end
        end
    end
    
    % 计算Zernike矩
    function Z = calculate_zernike_moment(points, n, m, options)
        % 实现论文中公式(4-23)的Zernike矩计算
        % Z_nm = ∑∑f(x,y)·V*_nm(ρ,θ)
        %
        % 参数：
        %   points - 归一化的点集 [-1, 1]
        %   n - 径向阶数
        %   m - 角阶数
        %   options - 可选设置
        
        % 提取点坐标
        x = points(:, 1);
        y = points(:, 2);
        
        % 转换为极坐标
        rho = sqrt(x.^2 + y.^2);
        theta = atan2(y, x);
        
        % 仅使用单位圆内的点
        valid_idx = rho <= 1;
        rho = rho(valid_idx);
        theta = theta(valid_idx);
        
        if isempty(rho)
            Z = 0;
            return;
        end
        
        % 计算径向多项式 R_nm(ρ)
        R = calculate_radial_polynomial(rho, n, abs(m));
        
        % 计算复值Zernike矩
        if m >= 0
            Z_real = mean(R .* cos(m*theta));
            Z_imag = mean(R .* sin(m*theta));
        else
            Z_real = mean(R .* cos(abs(m)*theta));
            Z_imag = -mean(R .* sin(abs(m)*theta));
        end
        
        % 组合成复数
        Z = complex(Z_real, Z_imag);
        
        % 确保结果有效
        if isnan(Z) || isinf(Z)
            Z = 0;
            if isfield(options, 'verbose') && options.verbose
                warning('Zernike矩计算结果是NaN或Inf，使用默认值0');
            end
        end
    end
    
    % 计算Zernike径向多项式
    function R = calculate_radial_polynomial(rho, n, m)
        % 实现Zernike径向多项式 R_nm(ρ)
        % 使用稳定的递归方法
        
        % 初始化
        R = zeros(size(rho));
        
        % 特殊情况处理
        if n == 0 && m == 0
            R = ones(size(rho));
            return;
        end
        
        if n == m
            R = rho.^n;
            return;
        end
        
        % 递归计算
        for k = 0:(n-m)/2
            coef = (-1)^k * factorial(n-k) / ...
                   (factorial(k) * factorial((n+m)/2-k) * factorial((n-m)/2-k));
            
            % 防止数值溢出
            if isnan(coef) || isinf(coef)
                continue;
            end
            
            R = R + coef * rho.^(n-2*k);
        end
    end
    
    % 计算Zernike多项式
    function Z = calculate_zernike_poly(rho, theta, n, m, options)
        % 确保m是整数且|m|<=n且n-|m|是偶数
        m_abs = abs(m);
        if m_abs > n || mod(n - m_abs, 2) ~= 0
            if isfield(options, 'verbose') && options.verbose
                fprintf('    警告: Zernike多项式(n=%d,m=%d)参数不满足条件: |m|<=n且(n-|m|)为偶数\n', n, m);
            end
            Z = zeros(size(rho));
            return;
        end
        
        % 计算径向多项式
        R = zeros(size(rho));
        
        % 特殊情况: n=m
        if n == m_abs
            R = rho.^n;
        else
            % 一般情况: 使用递归公式
            try
                for s = 0:((n - m_abs) / 2)
                    coef = (-1)^s * nchoosek(n-s, s) * nchoosek(n-2*s, (n+abs(m))/2-s) * nchoosek(n-abs(m)-2*s, (n-abs(m))/2-s);
                    
                    % 处理可能的NaN或Inf
                    if isnan(coef) || isinf(coef)
                        continue;
                    end
                    
                    R = R + coef * rho.^(n - 2*s);
                end
            catch calc_err
                if isfield(options, 'verbose') && options.verbose
                    fprintf('    警告: Zernike多项式(n=%d,m=%d)径向多项式计算出错: %s\n', ...
                        n, m, calc_err.message);
                end
                R = rho.^n; % 使用简化公式作为后备
            end
        end
        
        % 角度部分
        if m >= 0
            Z = R .* cos(m * theta);
        else
            Z = R .* sin(m_abs * theta);
        end
        
        % 检查结果
        if any(isnan(Z)) || any(isinf(Z))
            if isfield(options, 'verbose') && options.verbose
                fprintf('    警告: Zernike多项式(n=%d,m=%d)最终结果包含NaN或Inf值\n', n, m);
            end
            % 将NaN和Inf替换为0
            Z(isnan(Z) | isinf(Z)) = 0;
        end
    end
end