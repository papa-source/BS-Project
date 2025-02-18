% sine_plot.m
% �������Һ���ͼ��

% ���� x �����ݵ㣨�� 0 �� 4�У�
x = 0:0.01:4*pi;

% ��������ֵ
y = sin(x);

% ����ͼ��
figure;
plot(x, y, 'b-', 'LineWidth', 2);

% ��ӱ�������ǩ
title('���Һ���ͼ��');
xlabel('x');
ylabel('sin(x)');

% �������
grid on;

% ����������
axis([0 4*pi -1.2 1.2]);

% ���ˮƽ�ο���
hold on;
plot([0 4*pi], [0 0], 'k--', 'LineWidth', 0.5);
hold off; 

