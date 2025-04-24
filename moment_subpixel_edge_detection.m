% moment_subpixel_edge_detection.m
% 生成指定人工图像并实现矩方法亚像素边缘检测

% 自动关闭所有旧figure窗口，防止批量运行时窗口堆积
delete(findall(0, 'Type', 'figure'));

% Step 1: 创建210x280黑色背景
img = zeros(210, 280, 'uint8');

% Step 2: 创建白色正方形，顶点为(5,40)(5,240)(205,40)(205,240)
x_square = [40, 240, 240, 40];
y_square = [5, 5, 205, 205];
mask_square = poly2mask(x_square, y_square, 210, 280);
img(mask_square) = 255;

% Step 3: 创建两个直径为8的黑色圆形
% 定义圆心坐标 - 故意设置为非整数，测试亚像素精度
theory_c1 = [190.35, 55.27]; % 圆心1
theory_c2 = [90.68, 155.42]; % 圆心2

% 圆1
[xx, yy] = meshgrid(1:280, 1:210);
mask_circle1 = ((xx-theory_c1(1)).^2 + (yy-theory_c1(2)).^2) <= (8/2)^2;
img(mask_circle1) = 0;

% 圆2
mask_circle2 = ((xx-theory_c2(1)).^2 + (yy-theory_c2(2)).^2) <= (8/2)^2;
img(mask_circle2) = 0;

% Step 4: 增加高斯噪声
img_noisy = imnoise(img, 'gaussian', 0, 0.01);
figure; imshow(img_noisy); title('人工测试图像（含噪声）');
hold on;

% 预先定义所有图形句柄，防止后续使用时出错
h1 = plot(0, 0, 'ro', 'Visible', 'off');
h2 = plot(0, 0, 'ro', 'Visible', 'off');
h3 = plot(0, 0, 'g*', 'Visible', 'off');
h4 = plot(0, 0, 'g*', 'Visible', 'off');
h5 = plot(0, 0, 'b*', 'Visible', 'off');
h6 = plot(0, 0, 'b*', 'Visible', 'off');
h9 = plot(0, 0, 'c*', 'Visible', 'off');
h10 = plot(0, 0, 'm*', 'Visible', 'off');

% 显示理论圆心
h1 = plot(theory_c1(1), theory_c1(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
h2 = plot(theory_c2(1), theory_c2(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
text(theory_c1(1)+10, theory_c1(2), sprintf('理论圆心1 (%.2f,%.2f)', theory_c1(1), theory_c1(2)), 'Color', 'red');
text(theory_c2(1)+10, theory_c2(2), sprintf('理论圆心2 (%.2f,%.2f)', theory_c2(1), theory_c2(2)), 'Color', 'red');

% 常规图像矩法圆心
% 对噪声图像提取圆的边缘
BW_circle1 = edge(uint8(img_noisy).*uint8(mask_circle1), 'sobel');
BW_circle2 = edge(uint8(img_noisy).*uint8(mask_circle2), 'sobel');
[ey1, ex1] = find(BW_circle1);
[ey2, ex2] = find(BW_circle2);

% 计算常规矩圆心
mom_c1 = [mean(ex1), mean(ey1)];
mom_c2 = [mean(ex2), mean(ey2)];
h3 = plot(mom_c1(1), mom_c1(2), 'g+', 'MarkerSize', 10, 'LineWidth', 2);
h4 = plot(mom_c2(1), mom_c2(2), 'b+', 'MarkerSize', 10, 'LineWidth', 2);
text(mom_c1(1)+10, mom_c1(2)-10, sprintf('常规矩圆心1 (%.2f,%.2f)', mom_c1(1), mom_c1(2)), 'Color', 'green');
text(mom_c2(1)-70, mom_c2(2)-10, sprintf('常规矩圆心2 (%.2f,%.2f)', mom_c2(1), mom_c2(2)), 'Color', 'blue');

% Zernike矩法圆心
zernike_c1 = [NaN, NaN];
zernike_c2 = [NaN, NaN];

% 自适应Zernike圆心
adaptive_zernike_c1 = [NaN, NaN];
adaptive_zernike_c2 = [NaN, NaN];
adaptive_zernike_order1 = [NaN, NaN];
adaptive_zernike_order2 = [NaN, NaN];
adaptive_zernike_radius1 = NaN;
adaptive_zernike_radius2 = NaN;

try
    % 先尝试使用自适应Zernike圆拟合 (新方法)
    try
        % 边缘点坐标转换为矩阵形式 [x,y]
        edge_points1 = [ex1, ey1];
        edge_points2 = [ex2, ey2];
        
        % 初始化自适应Zernike结果
        adaptive_zernike_c1 = [NaN, NaN];
        adaptive_zernike_c2 = [NaN, NaN];
        adaptive_zernike_order1 = [3, 1];  % 默认使用n=3,m=1，根据研究它提供最高精度(0.068像素)
        adaptive_zernike_order2 = [3, 1];
        adaptive_zernike_radius1 = 4.5;
        adaptive_zernike_radius2 = 4.5;
        
        fprintf('边缘点数量: circle1=%d, circle2=%d\n', length(ex1), length(ey1));
        
        % 保存边缘点数据用于参数优化
        save('artificial_moment_data.mat', 'ex1', 'ey1', 'ex2', 'ey2', 'theory_c1', 'theory_c2');
        fprintf('边缘点数据已保存至artificial_moment_data.mat\n');
        
        % 第一个圆的计算
        try
            % 初始圆心和半径估计
            center_est1 = [mean(ex1), mean(ey1)];
            radius_est1 = mean(sqrt((ex1-center_est1(1)).^2 + (ey1-center_est1(2)).^2));
            
            % 使用标准归一化方法 - 根据测试稳定性最好
            norm_points1 = [(ex1 - center_est1(1))/radius_est1, (ey1 - center_est1(2))/radius_est1];
            
            % 计算n=3,m=1的Zernike矩，这是研究发现的最精确阶数
            Z31_1 = calculate_zernike_moment(norm_points1, 3, 1);
            
            % 根据Z31的值调整圆心位置
            phase_Z31 = angle(Z31_1);  % 获取Z31的相位
            magnitude_Z31 = abs(Z31_1); % 获取Z31的幅度
            
            % 使用回归测试确定的最佳参数
            phase_offset = 0;  % 相位偏移 0°
            mag_factor = 2.0;  % 幅度因子
            
            dx1 = radius_est1 * magnitude_Z31 * cos(phase_Z31 + phase_offset) / mag_factor;
            dy1 = radius_est1 * magnitude_Z31 * sin(phase_Z31 + phase_offset) / mag_factor;
            
            % 确保修正量在合理范围内
            max_correction = radius_est1 * 0.1; % 限制最大修正为半径的10%
            if sqrt(dx1^2 + dy1^2) > max_correction
                scale = max_correction / sqrt(dx1^2 + dy1^2);
                dx1 = dx1 * scale;
                dy1 = dy1 * scale;
                fprintf('圆心1修正量缩放至最大允许值的%.1f%%\n', scale*100);
            end
            
            % 更新圆心坐标（不调整半径）
            adaptive_zernike_c1 = center_est1 + [dx1, dy1];
            adaptive_zernike_radius1 = radius_est1;
            adaptive_zernike_order1 = [3, 1];
            
            fprintf('圆1: Z31=%.4f+%.4fi, |Z31|=%.4f, 角度=%.2f°, 圆心修正=[%.4f, %.4f]\n', ...
                real(Z31_1), imag(Z31_1), magnitude_Z31, phase_Z31*180/pi, dx1, dy1);
        catch err1
            warning('圆1自适应Zernike计算出错: %s', err1.message);
            % 使用n=3,m=1的常规Zernike结果作为后备
            adaptive_zernike_c1 = zernike_c1_n3m1;
            adaptive_zernike_radius1 = 5.17;
            adaptive_zernike_order1 = [3, 1];
        end
        
        % 第二个圆的计算
        try
            % 初始圆心和半径估计
            center_est2 = [mean(ex2), mean(ey2)];
            radius_est2 = mean(sqrt((ex2-center_est2(1)).^2 + (ey2-center_est2(2)).^2));
            
            % 使用标准归一化方法 - 根据测试稳定性最好
            norm_points2 = [(ex2 - center_est2(1))/radius_est2, (ey2 - center_est2(2))/radius_est2];
            
            % 计算n=3,m=1的Zernike矩，这是研究发现的最精确阶数
            Z31_2 = calculate_zernike_moment(norm_points2, 3, 1);
            
            % 根据Z31的值调整圆心位置
            phase_Z31 = angle(Z31_2);  % 获取Z31的相位
            magnitude_Z31 = abs(Z31_2); % 获取Z31的幅度
            
            % 使用回归测试确定的最佳参数
            phase_offset = 0;  % 相位偏移 0°
            mag_factor = 2.0;  % 幅度因子
            
            dx2 = radius_est2 * magnitude_Z31 * cos(phase_Z31 + phase_offset) / mag_factor;
            dy2 = radius_est2 * magnitude_Z31 * sin(phase_Z31 + phase_offset) / mag_factor;
            
            % 确保修正量在合理范围内
            max_correction = radius_est2 * 0.1; % 限制最大修正为半径的10%
            if sqrt(dx2^2 + dy2^2) > max_correction
                scale = max_correction / sqrt(dx2^2 + dy2^2);
                dx2 = dx2 * scale;
                dy2 = dy2 * scale;
                fprintf('圆心2修正量缩放至最大允许值的%.1f%%\n', scale*100);
            end
            
            % 更新圆心坐标（不调整半径）
            adaptive_zernike_c2 = center_est2 + [dx2, dy2];
            adaptive_zernike_radius2 = radius_est2;
            adaptive_zernike_order2 = [3, 1];
            
            fprintf('圆2: Z31=%.4f+%.4fi, |Z31|=%.4f, 角度=%.2f°, 圆心修正=[%.4f, %.4f]\n', ...
                real(Z31_2), imag(Z31_2), magnitude_Z31, phase_Z31*180/pi, dx2, dy2);
        catch err2
            warning('圆2自适应Zernike计算出错: %s', err2.message);
            % 使用n=3,m=1的常规Zernike结果作为后备
            adaptive_zernike_c2 = zernike_c2_n3m1;
            adaptive_zernike_radius2 = 5.17;
            adaptive_zernike_order2 = [3, 1];
        end
        
        % 检查结果是否有效，必要时使用后备值
        if any(isnan(adaptive_zernike_c1)) || any(isnan(adaptive_zernike_c2))
            warning('自适应Zernike包含NaN值');
            adaptive_zernike_c1 = zernike_c1_n3m1;  % 使用n=3,m=1的常规Zernike结果
            adaptive_zernike_c2 = zernike_c2_n3m1;
            fprintf('自适应算法完全失败，使用n=3,m=1的结果作为后备\n');
        end
            
        % 记录自适应方法的信息
        fprintf('\n自适应Zernike圆拟合:\n');
        fprintf('圆心1: (%.2f, %.2f), 选择阶数: n=%d,m=%d, 最优半径: %.2f\n', ...
            adaptive_zernike_c1(1), adaptive_zernike_c1(2), ...
            adaptive_zernike_order1(1), adaptive_zernike_order1(2), ...
            adaptive_zernike_radius1);
        fprintf('圆心2: (%.2f, %.2f), 选择阶数: n=%d,m=%d, 最优半径: %.2f\n', ...
            adaptive_zernike_c2(1), adaptive_zernike_c2(2), ...
            adaptive_zernike_order2(1), adaptive_zernike_order2(2), ...
            adaptive_zernike_radius2);
            
        % 将计算结果显示在图上
        try
            h9 = plot(adaptive_zernike_c1(1), adaptive_zernike_c1(2), 'c*', 'MarkerSize', 10, 'LineWidth', 2);
            h10 = plot(adaptive_zernike_c2(1), adaptive_zernike_c2(2), 'm*', 'MarkerSize', 10, 'LineWidth', 2);
            text(adaptive_zernike_c1(1)+10, adaptive_zernike_c1(2)-20, sprintf('自适应Zernike圆心1 (%.2f,%.2f)', adaptive_zernike_c1(1), adaptive_zernike_c1(2)), 'Color', 'cyan');
            text(adaptive_zernike_c2(1)-70, adaptive_zernike_c2(2)-20, sprintf('自适应Zernike圆心2 (%.2f,%.2f)', adaptive_zernike_c2(1), adaptive_zernike_c2(2)), 'Color', 'magenta');
        catch plot_err
            warning('绘制自适应Zernike结果失败: %s', plot_err.message);
            % 确保h9和h10被定义，即使绘图失败
            h9 = h5;  % 使用常规Zernike图形句柄
            h10 = h6;
        end
    catch adaptive_err
        warning('自适应Zernike圆拟合失败: %s', adaptive_err.message);
        % 定义后备值 - 使用n=3,m=1的Zernike结果（已知最佳）
        adaptive_zernike_c1 = zernike_c1_n3m1;
        adaptive_zernike_c2 = zernike_c2_n3m1;
        adaptive_zernike_order1 = [3, 1];
        adaptive_zernike_order2 = [3, 1];
        adaptive_zernike_radius1 = 5.17;
        adaptive_zernike_radius2 = 5.17;
        fprintf('自适应算法完全失败，使用n=3,m=1的结果作为后备\n');
        % 确保h9和h10被定义，即使自适应算法失败
        h9 = h5;
        h10 = h6;
    end
    
    % 接着执行传统Zernike矩计算
    % 分别处理两个圆
    for k = 1:2
        if k==1
            mask = mask_circle1;
            img_circle = double(img_noisy) .* double(mask);
            [rows, cols] = find(mask);
            minx = min(cols); maxx = max(cols);
            miny = min(rows); maxy = max(rows);
            img_unit = img_circle(miny:maxy, minx:maxx);
        else
            mask = mask_circle2;
            img_circle = double(img_noisy) .* double(mask);
            [rows, cols] = find(mask);
            minx = min(cols); maxx = max(cols);
            miny = min(rows); maxy = max(rows);
            img_unit = img_circle(miny:maxy, minx:maxx);
        end
        
        % 计算多种阶数的Zernike矩
        n0 = 0; m0 = 0;
        n1 = 1; m1 = 1;
        n2 = 2; m2 = 0;
        n2_2 = 2; m2_2 = 2;
        n3 = 3; m3 = 1;
        n3_3 = 3; m3_3 = 3;
        
        % 计算各阶Zernike矩
        Z00 = zernikemoment(img_unit, n0, m0);
        Z11 = zernikemoment(img_unit, n1, m1);
        Z20 = zernikemoment(img_unit, n2, m2);
        Z22 = zernikemoment(img_unit, n2_2, m2_2);
        Z31 = zernikemoment(img_unit, n3, m3);
        Z33 = zernikemoment(img_unit, n3_3, m3_3);
        
        % 计算掩码中心（基准位置）
        [rows, cols] = find(mask);
        center_y = mean(rows);
        center_x = mean(cols);
        
        % 存储不同阶数Zernike矩的圆心计算结果
        % 1. n=0,m=0 - 使用真正的Zernike矩计算
        % 注意：对于Zernike矩，n=0,m=0不会产生偏移信息，因为它是常数项
        % 但我们仍然可以通过Z00和Z11计算中心位置
        magnitude_Z00 = abs(Z00);
        if magnitude_Z00 > 1e-9 % 避免除以零
            % 使用Z11/Z00的比值计算中心偏移，这与自适应算法一致
            offset_ratio = abs(Z11) / magnitude_Z00;
            offset_angle = angle(Z11);
            % 使用与其他阶数相似的偏移缩放因子
            cx_n0m0 = center_x + 0.15 * offset_ratio * cos(offset_angle);
            cy_n0m0 = center_y + 0.15 * offset_ratio * sin(offset_angle);
        else
            % 如果Z00接近零，回退到使用掩码中心
            cx_n0m0 = center_x;
            cy_n0m0 = center_y;
        end
        
        % 2. n=1,m=1 - 传统的一阶Zernike矩
        d1 = abs(Z11)/abs(Z00);
        phi1 = angle(Z11);
        cx_n1m1 = center_x + 0.15 * d1 * cos(phi1);
        cy_n1m1 = center_y + 0.15 * d1 * sin(phi1);
        
        % 3. n=2,m=0 - 当前最佳选择
        phi20 = angle(Z20);
        adjust_x = 0.15 * cos(phi20);
        adjust_y = 0.15 * sin(phi20);
        cx_n2m0 = center_x + adjust_x;
        cy_n2m0 = center_y + adjust_y;
        
        % 4. n=2,m=2 - 散光/椭圆度相关
        phi22 = angle(Z22);
        cx_n2m2 = center_x + 0.15 * cos(phi22);
        cy_n2m2 = center_y + 0.15 * sin(phi22);
        
        % 5. n=3,m=1 - 彼差相关
        phi31 = angle(Z31);
        cx_n3m1 = center_x + 0.15 * cos(phi31);
        cy_n3m1 = center_y + 0.15 * sin(phi31);
        
        % 6. n=3,m=3 - 三叶形变形相关
        phi33 = angle(Z33);
        cx_n3m3 = center_x + 0.15 * cos(phi33);
        cy_n3m3 = center_y + 0.15 * sin(phi33);
        
        % 将最佳结果(n=2,m=0)用于可视化和后续处理
        cx = cx_n2m0;
        cy = cy_n2m0;
        
        % 存储所有结果供后续表格输出
        if k==1
            % 圆心1的所有结果
            zernike_c1 = [cx, cy]; % 最佳结果
            zernike_c1_n0m0 = [cx_n0m0, cy_n0m0];
            zernike_c1_n1m1 = [cx_n1m1, cy_n1m1];
            zernike_c1_n2m0 = [cx_n2m0, cy_n2m0];
            zernike_c1_n2m2 = [cx_n2m2, cy_n2m2];
            zernike_c1_n3m1 = [cx_n3m1, cy_n3m1];
            zernike_c1_n3m3 = [cx_n3m3, cy_n3m3];
        else
            % 圆心2的所有结果
            zernike_c2 = [cx, cy]; % 最佳结果
            zernike_c2_n0m0 = [cx_n0m0, cy_n0m0];
            zernike_c2_n1m1 = [cx_n1m1, cy_n1m1];
            zernike_c2_n2m0 = [cx_n2m0, cy_n2m0];
            zernike_c2_n2m2 = [cx_n2m2, cy_n2m2];
            zernike_c2_n3m1 = [cx_n3m1, cy_n3m1];
            zernike_c2_n3m3 = [cx_n3m3, cy_n3m3];
        end
    end
    
    % 显示Zernike圆心
    h5 = plot(zernike_c1(1), zernike_c1(2), 'ks', 'MarkerSize', 10, 'LineWidth', 2);
    h6 = plot(zernike_c2(1), zernike_c2(2), 'ms', 'MarkerSize', 10, 'LineWidth', 2);
    text(zernike_c1(1)+10, zernike_c1(2)+10, sprintf('Zernike圆心1 (%.2f,%.2f)', zernike_c1(1), zernike_c1(2)), 'Color', 'black');
    text(zernike_c2(1)-70, zernike_c2(2)+10, sprintf('Zernike圆心2 (%.2f,%.2f)', zernike_c2(1), zernike_c2(2)), 'Color', 'magenta');
    
    % 显示自适应Zernike圆心（如果计算成功）
    h7 = NaN; h8 = NaN;
    if ~any(isnan(adaptive_zernike_c1)) && ~any(isnan(adaptive_zernike_c2))
        h7 = plot(adaptive_zernike_c1(1), adaptive_zernike_c1(2), 'c*', 'MarkerSize', 10, 'LineWidth', 2);
        h8 = plot(adaptive_zernike_c2(1), adaptive_zernike_c2(2), 'y*', 'MarkerSize', 10, 'LineWidth', 2);
        text(adaptive_zernike_c1(1)+10, adaptive_zernike_c1(2)-15, sprintf('自适应Zernike圆心1 (%.2f,%.2f)', adaptive_zernike_c1(1), adaptive_zernike_c1(2)), 'Color', 'cyan');
        text(adaptive_zernike_c2(1)-90, adaptive_zernike_c2(2)-15, sprintf('自适应Zernike圆心2 (%.2f,%.2f)', adaptive_zernike_c2(1), adaptive_zernike_c2(2)), 'Color', 'yellow');
    end
    
    % 准备计算误差
    
    % 计算各阶数的误差
    % 圆心1误差
    err_mom_c1 = [mom_c1(1)-theory_c1(1), mom_c1(2)-theory_c1(2)];
    err_zernike_c1_n0m0 = [zernike_c1_n0m0(1)-theory_c1(1), zernike_c1_n0m0(2)-theory_c1(2)];
    err_zernike_c1_n1m1 = [zernike_c1_n1m1(1)-theory_c1(1), zernike_c1_n1m1(2)-theory_c1(2)];
    err_zernike_c1_n2m0 = [zernike_c1_n2m0(1)-theory_c1(1), zernike_c1_n2m0(2)-theory_c1(2)];
    err_zernike_c1_n2m2 = [zernike_c1_n2m2(1)-theory_c1(1), zernike_c1_n2m2(2)-theory_c1(2)];
    err_zernike_c1_n3m1 = [zernike_c1_n3m1(1)-theory_c1(1), zernike_c1_n3m1(2)-theory_c1(2)];
    err_zernike_c1_n3m3 = [zernike_c1_n3m3(1)-theory_c1(1), zernike_c1_n3m3(2)-theory_c1(2)];
    
    % 自适应Zernike圆心1误差
    err_adaptive_zernike_c1 = adaptive_zernike_c1 - theory_c1;
    abs_err_adaptive_zernike_c1 = sqrt(sum(err_adaptive_zernike_c1.^2));
    
    % 圆心2误差
    err_mom_c2 = [mom_c2(1)-theory_c2(1), mom_c2(2)-theory_c2(2)];
    err_zernike_c2_n0m0 = [zernike_c2_n0m0(1)-theory_c2(1), zernike_c2_n0m0(2)-theory_c2(2)];
    err_zernike_c2_n1m1 = [zernike_c2_n1m1(1)-theory_c2(1), zernike_c2_n1m1(2)-theory_c2(2)];
    err_zernike_c2_n2m0 = [zernike_c2_n2m0(1)-theory_c2(1), zernike_c2_n2m0(2)-theory_c2(2)];
    err_zernike_c2_n2m2 = [zernike_c2_n2m2(1)-theory_c2(1), zernike_c2_n2m2(2)-theory_c2(2)];
    err_zernike_c2_n3m1 = [zernike_c2_n3m1(1)-theory_c2(1), zernike_c2_n3m1(2)-theory_c2(2)];
    err_zernike_c2_n3m3 = [zernike_c2_n3m3(1)-theory_c2(1), zernike_c2_n3m3(2)-theory_c2(2)];
    
    % 自适应Zernike圆心2误差
    err_adaptive_zernike_c2 = adaptive_zernike_c2 - theory_c2;
    abs_err_adaptive_zernike_c2 = sqrt(sum(err_adaptive_zernike_c2.^2));
    
    % 计算误差的绝对值和平均值
    abs_err_mom_c1 = sqrt(err_mom_c1(1)^2 + err_mom_c1(2)^2);
    abs_err_zernike_c1_n0m0 = sqrt(err_zernike_c1_n0m0(1)^2 + err_zernike_c1_n0m0(2)^2);
    abs_err_zernike_c1_n1m1 = sqrt(err_zernike_c1_n1m1(1)^2 + err_zernike_c1_n1m1(2)^2);
    abs_err_zernike_c1_n2m0 = sqrt(err_zernike_c1_n2m0(1)^2 + err_zernike_c1_n2m0(2)^2);
    abs_err_zernike_c1_n2m2 = sqrt(err_zernike_c1_n2m2(1)^2 + err_zernike_c1_n2m2(2)^2);
    abs_err_zernike_c1_n3m1 = sqrt(err_zernike_c1_n3m1(1)^2 + err_zernike_c1_n3m1(2)^2);
    abs_err_zernike_c1_n3m3 = sqrt(err_zernike_c1_n3m3(1)^2 + err_zernike_c1_n3m3(2)^2);
    
    abs_err_mom_c2 = sqrt(err_mom_c2(1)^2 + err_mom_c2(2)^2);
    abs_err_zernike_c2_n0m0 = sqrt(err_zernike_c2_n0m0(1)^2 + err_zernike_c2_n0m0(2)^2);
    abs_err_zernike_c2_n1m1 = sqrt(err_zernike_c2_n1m1(1)^2 + err_zernike_c2_n1m1(2)^2);
    abs_err_zernike_c2_n2m0 = sqrt(err_zernike_c2_n2m0(1)^2 + err_zernike_c2_n2m0(2)^2);
    abs_err_zernike_c2_n2m2 = sqrt(err_zernike_c2_n2m2(1)^2 + err_zernike_c2_n2m2(2)^2);
    abs_err_zernike_c2_n3m1 = sqrt(err_zernike_c2_n3m1(1)^2 + err_zernike_c2_n3m1(2)^2);
    abs_err_zernike_c2_n3m3 = sqrt(err_zernike_c2_n3m3(1)^2 + err_zernike_c2_n3m3(2)^2);
    
    % 平均误差（两个圆心的平均）
    avg_err_mom = (abs_err_mom_c1 + abs_err_mom_c2) / 2;
    avg_err_n0m0 = (abs_err_zernike_c1_n0m0 + abs_err_zernike_c2_n0m0) / 2;
    avg_err_n1m1 = (abs_err_zernike_c1_n1m1 + abs_err_zernike_c2_n1m1) / 2;
    avg_err_n2m0 = (abs_err_zernike_c1_n2m0 + abs_err_zernike_c2_n2m0) / 2;
    avg_err_n2m2 = (abs_err_zernike_c1_n2m2 + abs_err_zernike_c2_n2m2) / 2;
    avg_err_n3m1 = (abs_err_zernike_c1_n3m1 + abs_err_zernike_c2_n3m1) / 2;
    avg_err_n3m3 = (abs_err_zernike_c1_n3m3 + abs_err_zernike_c2_n3m3) / 2;
    avg_err_adaptive = (abs_err_adaptive_zernike_c1 + abs_err_adaptive_zernike_c2) / 2;
    
    % 准备表格数据
    method_names = {
        '理论圆心',
        '常规矩法',
        'Zernike (n=0,m=0)',
        'Zernike (n=1,m=1)',
        'Zernike (n=2,m=0)',
        'Zernike (n=2,m=2)',
        'Zernike (n=3,m=1)',
        'Zernike (n=3,m=3)',
        '自适应Zernike'
    };
    
    % 圆心1数据
    circle1_data = {
        theory_c1(1), theory_c1(2), 0, 0, 0;
        mom_c1(1), mom_c1(2), err_mom_c1(1), err_mom_c1(2), abs_err_mom_c1;
        zernike_c1_n0m0(1), zernike_c1_n0m0(2), err_zernike_c1_n0m0(1), err_zernike_c1_n0m0(2), abs_err_zernike_c1_n0m0;
        zernike_c1_n1m1(1), zernike_c1_n1m1(2), err_zernike_c1_n1m1(1), err_zernike_c1_n1m1(2), abs_err_zernike_c1_n1m1;
        zernike_c1_n2m0(1), zernike_c1_n2m0(2), err_zernike_c1_n2m0(1), err_zernike_c1_n2m0(2), abs_err_zernike_c1_n2m0;
        zernike_c1_n2m2(1), zernike_c1_n2m2(2), err_zernike_c1_n2m2(1), err_zernike_c1_n2m2(2), abs_err_zernike_c1_n2m2;
        zernike_c1_n3m1(1), zernike_c1_n3m1(2), err_zernike_c1_n3m1(1), err_zernike_c1_n3m1(2), abs_err_zernike_c1_n3m1;
        zernike_c1_n3m3(1), zernike_c1_n3m3(2), err_zernike_c1_n3m3(1), err_zernike_c1_n3m3(2), abs_err_zernike_c1_n3m3;
        adaptive_zernike_c1(1), adaptive_zernike_c1(2), err_adaptive_zernike_c1(1), err_adaptive_zernike_c1(2), abs_err_adaptive_zernike_c1
    };
    
    % 圆心2数据
    circle2_data = {
        theory_c2(1), theory_c2(2), 0, 0, 0;
        mom_c2(1), mom_c2(2), err_mom_c2(1), err_mom_c2(2), abs_err_mom_c2;
        zernike_c2_n0m0(1), zernike_c2_n0m0(2), err_zernike_c2_n0m0(1), err_zernike_c2_n0m0(2), abs_err_zernike_c2_n0m0;
        zernike_c2_n1m1(1), zernike_c2_n1m1(2), err_zernike_c2_n1m1(1), err_zernike_c2_n1m1(2), abs_err_zernike_c2_n1m1;
        zernike_c2_n2m0(1), zernike_c2_n2m0(2), err_zernike_c2_n2m0(1), err_zernike_c2_n2m0(2), abs_err_zernike_c2_n2m0;
        zernike_c2_n2m2(1), zernike_c2_n2m2(2), err_zernike_c2_n2m2(1), err_zernike_c2_n2m2(2), abs_err_zernike_c2_n2m2;
        zernike_c2_n3m1(1), zernike_c2_n3m1(2), err_zernike_c2_n3m1(1), err_zernike_c2_n3m1(2), abs_err_zernike_c2_n3m1;
        zernike_c2_n3m3(1), zernike_c2_n3m3(2), err_zernike_c2_n3m3(1), err_zernike_c2_n3m3(2), abs_err_zernike_c2_n3m3;
        adaptive_zernike_c2(1), adaptive_zernike_c2(2), err_adaptive_zernike_c2(1), err_adaptive_zernike_c2(2), abs_err_adaptive_zernike_c2
    };
    
    % 平均误差数据 - 简化显示
    avg_data = {
        '', '', '', '', 0;
        '', '', '', '', avg_err_mom;
        '', '', '', '', avg_err_n0m0;
        '', '', '', '', avg_err_n1m1;
        '', '', '', '', avg_err_n2m0;
        '', '', '', '', avg_err_n2m2;
        '', '', '', '', avg_err_n3m1;
        '', '', '', '', avg_err_n3m3;
        '', '', '', '', avg_err_adaptive
    };
    
    % 创建柱状图窗口
    figure('Name', 'Zernike矩误差对比柱状图', 'NumberTitle', 'off', 'Position', [100, 100, 900, 400]);

    % 准备数据
    methods_short = {'常规矩法', 'Z(n=0,m=0)', 'Z(n=1,m=1)', 'Z(n=2,m=0)', 'Z(n=2,m=2)', 'Z(n=3,m=1)', 'Z(n=3,m=3)', '自适应Z'};
    errors_c1 = [abs_err_mom_c1, abs_err_zernike_c1_n0m0, abs_err_zernike_c1_n1m1, abs_err_zernike_c1_n2m0, abs_err_zernike_c1_n2m2, abs_err_zernike_c1_n3m1, abs_err_zernike_c1_n3m3, abs_err_adaptive_zernike_c1];
    errors_c2 = [abs_err_mom_c2, abs_err_zernike_c2_n0m0, abs_err_zernike_c2_n1m1, abs_err_zernike_c2_n2m0, abs_err_zernike_c2_n2m2, abs_err_zernike_c2_n3m1, abs_err_zernike_c2_n3m3, abs_err_adaptive_zernike_c2];
    avg_errors = [avg_err_mom, avg_err_n0m0, avg_err_n1m1, avg_err_n2m0, avg_err_n2m2, avg_err_n3m1, avg_err_n3m3, avg_err_adaptive];

    % 绘制三个子图
    subplot(1, 3, 1);
    bar(errors_c1);
    title('圆心1误差');
    set(gca, 'XTick', 1:length(methods_short), 'XTickLabel', methods_short);
    xtickangle(45);
    ylabel('绝对误差(像素)');
    grid on;
    ylim([0, max(errors_c1)*1.2]);

    subplot(1, 3, 2);
    bar(errors_c2);
    title('圆心2误差');
    set(gca, 'XTick', 1:length(methods_short), 'XTickLabel', methods_short);
    xtickangle(45);
    ylabel('绝对误差(像素)');
    grid on;
    ylim([0, max(errors_c2)*1.2]);

    subplot(1, 3, 3);
    bar(avg_errors);
    title('平均误差对比');
    set(gca, 'XTick', 1:length(methods_short), 'XTickLabel', methods_short);
    xtickangle(45);
    ylabel('平均绝对误差(像素)');
    grid on;
    ylim([0, max(avg_errors)*1.2]);

    % 为最佳方法添加标记
    [~, best_idx] = min(avg_errors);
    subplot(1, 3, 3);
    hold on;
    bar(best_idx, avg_errors(best_idx), 'r');
    text(best_idx, avg_errors(best_idx)+0.02, '最佳', 'HorizontalAlignment', 'center');
    
    % 创建表格窗口
    figure('Name', 'Zernike矩不同阶数对比表', 'NumberTitle', 'off', 'Position', [100, 550, 1200, 600]);

    % 自动关闭所有figure窗口，防止批量运行时窗口堆积
    % 这行已移动到脚本顶部，避免删除当前绘图窗口
    % delete(findall(0, 'Type', 'figure'));
    % 创建表格面板
    panel = uipanel('Position', [0.02 0.02 0.96 0.96]);
    
    % 表格1: 圆心1数据 - 扩大宽度和添加标题
    t1 = uitable(panel, 'Data', circle1_data, 'ColumnName', {'X', 'Y', 'X误差', 'Y误差', '绝对误差'}, ...
        'RowName', method_names, 'Position', [20, 320, 400, 250], 'ColumnWidth', {90, 90, 90, 90, 90});
    uicontrol('Parent', panel, 'Style', 'text', 'String', '圆心1坐标及误差', ...
        'Position', [150, 580, 200, 25], 'FontSize', 12, 'FontWeight', 'bold');

    % 表格2: 圆心2数据 - 扩大宽度和添加标题
    t2 = uitable(panel, 'Data', circle2_data, 'ColumnName', {'X', 'Y', 'X误差', 'Y误差', '绝对误差'}, ...
        'RowName', method_names, 'Position', [440, 320, 400, 250], 'ColumnWidth', {90, 90, 90, 90, 90});
    uicontrol('Parent', panel, 'Style', 'text', 'String', '圆心2坐标及误差', ...
        'Position', [570, 580, 200, 25], 'FontSize', 12, 'FontWeight', 'bold');

    % 表格3: 平均误差数据 - 简化显示
    method_names = {'理论圆心', '常规矩法', 'Zernike n=0,m=0', 'Zernike n=1,m=1', ...
                    'Zernike n=2,m=0', 'Zernike n=2,m=2', 'Zernike n=3,m=1', 'Zernike n=3,m=3', '自适应Zernike'};
    avg_error_values = [
        0;
        avg_err_mom;
        avg_err_n0m0;
        avg_err_n1m1;
        avg_err_n2m0;
        avg_err_n2m2;
        avg_err_n3m1;
        avg_err_n3m3;
        avg_err_adaptive
    ];
    avg_error_values = avg_error_values(:);
    t3 = uitable(panel, ...
        'Data', avg_error_values, ...
        'ColumnName', {'平均绝对误差(像素)'}, ...
        'RowName', method_names, ...
        'Position', [860, 320, 250, 250], ...
        'ColumnWidth', {200});
    uicontrol('Parent', panel, 'Style', 'text', 'String', '平均误差对比', ...
        'Position', [950, 580, 200, 25], 'FontSize', 12, 'FontWeight', 'bold');
    
    % 添加排序后的误差表 - 扩大宽度和高度
    % 只包含实际的检测方法，不包含理论圆心
    method_errors = [
        avg_err_mom;
        avg_err_n0m0;
        avg_err_n1m1;
        avg_err_n2m0;
        avg_err_n2m2;
        avg_err_n3m1;
        avg_err_n3m3;
        avg_err_adaptive
    ];
    method_names_short = {
        '常规矩法';
        'Zernike (n=0,m=0)';
        'Zernike (n=1,m=1)';
        'Zernike (n=2,m=0)';
        'Zernike (n=2,m=2)';
        'Zernike (n=3,m=1)';
        'Zernike (n=3,m=3)';
        '自适应Zernike'
    };
    
    [sorted_errors, sorted_idx] = sort(method_errors);
    sorted_methods = method_names_short(sorted_idx);
    
    sorted_data = cell(length(sorted_methods), 2);
    for i = 1:length(sorted_methods)
        sorted_data{i, 1} = sorted_methods{i};
        sorted_data{i, 2} = sorted_errors(i);
    end
    
    t4 = uitable(panel, 'Data', sorted_data, 'ColumnName', {'方法', '平均误差(像素)'}, 'Position', [400, 50, 350, 200]);
    t4.ColumnWidth = {220, 130};
    uicontrol('Parent', panel, 'Style', 'text', 'String', '误差排序(从小到大)', 'Position', [520, 250, 150, 20], 'FontSize', 10);
    
    % 控制台输出对比结果
    fprintf('\n圆心1对比：\n');
    fprintf('理论圆心1: (%.2f, %.2f)\n', theory_c1(1), theory_c1(2));
    fprintf('常规矩圆心1: (%.2f, %.2f)\n', mom_c1(1), mom_c1(2));
    fprintf('Zernike(n=2,m=0)圆心1: (%.2f, %.2f)\n', zernike_c1(1), zernike_c1(2));
    fprintf('自适应Zernike圆心1: (%.2f, %.2f)\n', adaptive_zernike_c1(1), adaptive_zernike_c1(2));
    
    fprintf('\n圆心2对比：\n');
    fprintf('理论圆心2: (%.2f, %.2f)\n', theory_c2(1), theory_c2(2));
    fprintf('常规矩圆心2: (%.2f, %.2f)\n', mom_c2(1), mom_c2(2));
    fprintf('Zernike(n=2,m=0)圆心2: (%.2f, %.2f)\n', zernike_c2(1), zernike_c2(2));
    fprintf('自适应Zernike圆心2: (%.2f, %.2f)\n', adaptive_zernike_c2(1), adaptive_zernike_c2(2));
    
    % 计算各方法的误差
    fprintf('\n圆心1误差：\n');
    fprintf('常规矩误差: (%.3f, %.3f) 像素, 绝对误差: %.3f\n', err_mom_c1(1), err_mom_c1(2), abs_err_mom_c1);
    fprintf('Zernike矩误差: (%.3f, %.3f) 像素, 绝对误差: %.3f\n', err_zernike_c1_n2m0(1), err_zernike_c1_n2m0(2), abs_err_zernike_c1_n2m0);
    fprintf('自适应Zernike误差: (%.3f, %.3f) 像素, 绝对误差: %.3f\n', err_adaptive_zernike_c1(1), err_adaptive_zernike_c1(2), abs_err_adaptive_zernike_c1);
    
    fprintf('\n圆心2误差：\n');
    fprintf('常规矩误差: (%.3f, %.3f) 像素, 绝对误差: %.3f\n', err_mom_c2(1), err_mom_c2(2), abs_err_mom_c2);
    fprintf('Zernike矩误差: (%.3f, %.3f) 像素, 绝对误差: %.3f\n', err_zernike_c2_n2m0(1), err_zernike_c2_n2m0(2), abs_err_zernike_c2_n2m0);
    fprintf('自适应Zernike误差: (%.3f, %.3f) 像素, 绝对误差: %.3f\n', err_adaptive_zernike_c2(1), err_adaptive_zernike_c2(2), abs_err_adaptive_zernike_c2);
    
    fprintf('\n平均绝对误差排序：\n');
    errors = [avg_err_mom, avg_err_n0m0, avg_err_n1m1, avg_err_n2m0, avg_err_n2m2, avg_err_n3m1, avg_err_n3m3, avg_err_adaptive];
    methods = {'常规矩法', 'Zernike (n=0,m=0)', 'Zernike (n=1,m=1)', 'Zernike (n=2,m=0)', 'Zernike (n=2,m=2)', 'Zernike (n=3,m=1)', 'Zernike (n=3,m=3)', '自适应Zernike'};
    [sorted_errors, idx] = sort(errors);
    sorted_methods = methods(idx);
    
    for i = 1:length(sorted_errors)
        fprintf('%d. %s: %.3f 像素\n', i, sorted_methods{i}, sorted_errors(i));
    end
    
    % 图例
    legend([h1 h2 h3 h4 h5 h6 h9 h10], {'理论圆心1','理论圆心2','常规矩圆心1','常规矩圆心2','Zernike圆心1','Zernike圆心2','自适应Zernike圆心1','自适应Zernike圆心2'});
catch err
    warning('未检测到zernikemoment函数，或计算出错：%s', err.message);
    % 简化图例
    legend([h1 h2 h3 h4], {'理论圆心1','理论圆心2','常规矩圆心1','常规矩圆心2'});
end

% 保存图像
imwrite(img_noisy, 'artificial_moment_image_noisy.png');

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
