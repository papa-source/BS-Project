% sine_plot.m
% 绘制正弦函数图像

% 创建 x 轴数据点（从 0 到 4π）
x = 0:0.01:4*pi;

% 计算正弦值
y = sin(x);

% 创建图像
figure;
plot(x, y, 'b-', 'LineWidth', 2);

% 添加标题和轴标签
title('正弦函数图像');
xlabel('x');
ylabel('sin(x)');

% 添加网格
grid on;

% 设置坐标轴
axis([0 4*pi -1.2 1.2]);

% 添加水平参考线
hold on;
plot([0 4*pi], [0 0], 'k--', 'LineWidth', 0.5);
hold off; 

