function display_results(original_img, edges, centers, radii, lines, corners)
% ��ʾ�����
figure('Name', '��ĳ�����������', 'NumberTitle', 'off');

% 1. ԭͼ
subplot(2,2,1);
imshow(original_img);
title('ԭʼͼ��', 'FontSize', 12);

% 2. ��Ե�����
subplot(2,2,2);
imshow(edges);
title('��Ե�����', 'FontSize', 12);

% 3. Բ�׼����
subplot(2,2,3);
imshow(original_img);
hold on;
if ~isempty(centers)
    % �Ȼ�Բ
    viscircles(centers, radii, 'EdgeColor', 'b', 'LineWidth', 1.5);
    
    % ��Բ�׽������򣨴����ϵ����£�
    [~, order] = sortrows(centers, [2 1]);  % �Ȱ�y�����ٰ�x����
    
    % ��עԲ�ĺͰ뾶��ʹ�ø������Ĳ���
    for i = 1:size(centers, 1)
        idx = order(i);
        % ��Բ��
        plot(centers(idx,1), centers(idx,2), 'b+', 'MarkerSize', 8, 'LineWidth', 1.5);
        
        % �����עλ�ã������ص���
        angle = 45;  % ��עλ�õĽǶ�
        offset = radii(idx) * 1.2;  % ��ע����Բ�ĵľ���
        text_x = centers(idx,1) + offset * cosd(angle);
        text_y = centers(idx,2) - offset * sind(angle);
        
        % ��ӱ�ע�ı�
        text(text_x, text_y, sprintf('%d\nR=%.1f', i, radii(idx)), ...
             'Color', 'blue', 'FontSize', 8, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'left', ...
             'BackgroundColor', [1 1 1 0.7]);  % ��͸����ɫ����
    end
end
title('Բ�׼����', 'FontSize', 12);
hold off;

% 4. �����߼����
subplot(2,2,4);
imshow(original_img);
hold on;

% ����������
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1), xy(:,2), 'LineWidth', 2, 'Color', [0 0.7 0]);
end

% ���ƽǵ�
if ~isempty(corners)
    plot(corners(:,1), corners(:,2), 'r*', 'MarkerSize', 8);
    % ��ע�ǵ���
    for i = 1:size(corners, 1)
        text(corners(i,1)+5, corners(i,2)+5, sprintf('%d', i), ...
             'Color', 'red', 'FontSize', 8, 'FontWeight', 'bold', ...
             'BackgroundColor', [1 1 1 0.7]);
    end
end

title('�����߼����', 'FontSize', 12);
hold off;

% ����ͼ����ʾ
set(gcf, 'Position', get(0, 'Screensize'));  % ȫ����ʾ
end