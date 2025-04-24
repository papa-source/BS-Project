% 批量运行 Zernike 矩检测主程序，统计最小误差分布并保存详细误差数据
N = 10; % 运行次数
% 定义所有检测方法的名称（不包含理论圆心）
% 只包含实际检测方法（7个），与err_matrix列数严格一致
method_names = {'常规矩法', 'Zernike n=0,m=0', 'Zernike n=1,m=1', 'Zernike n=2,m=0', 'Zernike n=2,m=2', 'Zernike n=3,m=1', 'Zernike n=3,m=3'};
winner_count = zeros(1, length(method_names));
err_matrix = zeros(N, length(method_names)); % 每次所有方法的误差
min_errs = zeros(N, 1);                    % 每次实验最小误差
min_idx = zeros(N, 1);                     % 每次最优方法索引

for i = 1:N
    run('moment_subpixel_edge_detection.m');
    % 确保误差数据和方法名称的顺序一致
    errs = [avg_err_mom, avg_err_n0m0, avg_err_n1m1, avg_err_n2m0, avg_err_n2m2, avg_err_n3m1, avg_err_n3m3];
    err_matrix(i, :) = errs;
    [min_errs(i), min_idx(i)] = min(errs);
    winner_count(min_idx(i)) = winner_count(min_idx(i)) + 1;
end

% 保存详细误差数据
save('zernike_batch_results.mat', 'err_matrix', 'min_errs', 'min_idx', 'winner_count', 'method_names');

% 可视化1：各方法最优次数
figure;
bar(winner_count);
set(gca, 'XTickLabel', method_names, 'XTickLabelRotation', 30);
ylabel('成为最优的次数');
title(['各阶Zernike矩成为最优的统计次数（', num2str(N), '次实验）']);
grid on;

% 可视化2：最小误差分布直方图
figure;
histogram(min_errs, 20);
xlabel('每次实验的最小误差值');
ylabel('出现次数');
title(['最小误差值分布直方图（', num2str(N), '次实验）']);
grid on;

% 可视化3：所有方法误差分布箱线图
figure;

% 强制重新设置method_names为正确的7个方法，避免与主程序中的全局变量冲突
box_labels = {'常规矩法', 'Zernike n=0,m=0', 'Zernike n=1,m=1', 'Zernike n=2,m=0', 'Zernike n=2,m=2', 'Zernike n=3,m=1', 'Zernike n=3,m=3'};

% 严格确保标签数量与数据列数匹配
if size(err_matrix, 2) ~= length(box_labels)
    fprintf('流程信息: 数据列数 %d, 标签数量 %d\n', size(err_matrix, 2), length(box_labels));
    % 确保数据和标签长度匹配
    box_labels = box_labels(1:size(err_matrix, 2));
    fprintf('已自动调整标签数量为 %d\n', length(box_labels));
end

% 使用正确数量的标签调用boxplot
boxplot(err_matrix, 'Labels', box_labels);
ylabel('平均绝对误差（像素）');
title(['各方法误差分布箱线图（', num2str(N), '次实验）']);
grid on;

% 关闭可能的空白图或多余窗口
figs = findall(0, 'Type', 'figure');
keep_figs = figs(end-2:end);
del_figs = setdiff(figs, keep_figs);
for f = del_figs'
    try 
        delete(f);
    catch
        % 忽略删除错误
    end
end
