function [centers, radii] = detect7(edges)
% Բ��ʶ�������
% ʹ�ñ�Ե���������Բ�׼������

% ʹ�ñ�Ե���������Բ���
[centers, radii] = imfindcircles(edges, [10 30], ...  % �����뾶��ΧΪ10-30����
    'ObjectPolarity', 'bright', ...
    'Sensitivity', 0.85, ...  % �������ж�
    'EdgeThreshold', 0.1, ...
    'Method', 'PhaseCode');

% ���û�м�⵽Բ������ʹ�ö�ֵ����
if isempty(centers)
    % ע�⣺������Ҫ���÷��ṩholes_mask�����������ں�����������Ӹò���
    % Ϊ�˱��ּ����ԣ��������Ƿ��ؿս��
    centers = [];
    radii = [];
end
end