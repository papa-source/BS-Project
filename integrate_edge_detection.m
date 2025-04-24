function edges = integrate_edge_detection(img, use_gui)
% INTEGRATE_EDGE_DETECTION 边缘检测集成函数
% 该函数集成多种边缘检测算法，可以通过GUI或直接调用不同算法
% 输入:
%   img - 输入图像 (灰度图)
%   use_gui - 是否使用GUI界面选择算法 (默认为true)
% 输出:
%   edges - 边缘检测结果

% 参数检查
if nargin < 2
    use_gui = true;
end

% 确保图像是灰度图
if size(img, 3) == 3
    img = rgb2gray(img);
end

% 方法1: 使用GUI界面选择算法
if use_gui
    % 打开GUI界面
    edge_detection_gui();
    
    % 将图像传递给GUI并等待处理
    assignin('base', 'currentImage', img);
    
    % 显示提示消息
    h = msgbox('请在GUI中选择边缘检测算法和参数，然后点击"应用到主程序"按钮。\n\n当完成设置后，点击此消息框的"OK"按钮继续。', '提示', 'modal');
    uiwait(h);
    
    % 检查是否生成了结果
    if evalin('base', 'exist(''edgeDetectionResult'', ''var'')')
        edgeResult = evalin('base', 'edgeDetectionResult');
        edges = edgeResult.image;
        
        % 显示所选算法信息
        disp(['已选择 ', edgeResult.method, ' 边缘检测算法']);
        
        % 输出参数信息
        switch edgeResult.method
            case 'Canny'
                disp(['参数: 低阈值=', num2str(edgeResult.params.low), ', 高阈值=', num2str(edgeResult.params.high), ', Sigma=', num2str(edgeResult.params.sigma)]);
            case {'Sobel', 'Prewitt', 'Roberts'}
                disp(['参数: 阈值=', num2str(edgeResult.params.threshold)]);
            case 'LoG'
                disp(['参数: 阈值=', num2str(edgeResult.params.threshold), ', Sigma=', num2str(edgeResult.params.sigma)]);
        end
    else
        % 如果用户未应用结果，则使用默认Canny算法
        disp('未检测到边缘检测结果，使用默认Canny算法');
        edges = edge(img, 'Canny');
    end
    
% 方法2: 自动选择算法 - 如不使用GUI界面
else
    % 弹出对话框让用户选择算法
    methods = {'Canny', 'Sobel', 'Prewitt', 'Roberts', 'LoG', '比较所有算法'};
    [selection, ok] = listdlg('ListString', methods, 'SelectionMode', 'single', ...
        'Name', '选择边缘检测算法', 'PromptString', '请选择一种边缘检测算法:', ...
        'ListSize', [200 150]);
    
    if ok
        switch selection
            case 1 % Canny
                % 使用自动确定的阈值
                [~, threshold] = edge(img, 'Canny');
                edges = edge(img, 'Canny', threshold);
                disp(['使用Canny算法，阈值: [', num2str(threshold(1)), ', ', num2str(threshold(2)), ']']);
                
            case 2 % Sobel
                threshold = graythresh(img) * 0.5; % 使用Otsu阈值的一半作为初始值
                edges = edge(img, 'Sobel', threshold);
                disp(['使用Sobel算法，阈值: ', num2str(threshold)]);
                
            case 3 % Prewitt
                threshold = graythresh(img) * 0.5;
                edges = edge(img, 'Prewitt', threshold);
                disp(['使用Prewitt算法，阈值: ', num2str(threshold)]);
                
            case 4 % Roberts
                threshold = graythresh(img) * 0.5;
                edges = edge(img, 'Roberts', threshold);
                disp(['使用Roberts算法，阈值: ', num2str(threshold)]);
                
            case 5 % LoG
                threshold = 0.003;
                sigma = 2.0;
                edges = edge(img, 'log', threshold, sigma);
                disp(['使用LoG算法，阈值: ', num2str(threshold), ', Sigma: ', num2str(sigma)]);
                
            case 6 % 比较所有算法
                % 调用比较函数
                [edges_canny, ~, edges_canny_final, edges_roberts, edges_sobel, edges_prewitt, edges_log] = compare_edge_detectors(img);
                
                % 弹出对话框让用户选择最终使用的算法
                results = {'优化后的Canny', '原始Canny', 'Roberts', 'Sobel', 'Prewitt', 'LoG'};
                [selection2, ok2] = listdlg('ListString', results, 'SelectionMode', 'single', ...
                    'Name', '选择最终结果', 'PromptString', '请选择要应用的最终边缘检测结果:', ...
                    'ListSize', [200 150]);
                
                if ok2
                    switch selection2
                        case 1
                            edges = edges_canny_final;
                            disp('使用优化后的Canny算法结果');
                        case 2
                            edges = edges_canny;
                            disp('使用原始Canny算法结果');
                        case 3
                            edges = edges_roberts;
                            disp('使用Roberts算法结果');
                        case 4
                            edges = edges_sobel;
                            disp('使用Sobel算法结果');
                        case 5
                            edges = edges_prewitt;
                            disp('使用Prewitt算法结果');
                        case 6
                            edges = edges_log;
                            disp('使用LoG算法结果');
                    end
                else
                    % 默认使用优化后的Canny
                    edges = edges_canny_final;
                    disp('用户取消选择，默认使用优化后的Canny算法结果');
                end
        end
    else
        % 用户取消，使用默认Canny
        edges = edge(img, 'Canny');
        disp('用户取消选择，使用默认Canny算法');
    end
end

end

% 从工作区加载图像的辅助函数 (供GUI使用)
function loadImageFromWS()
    if evalin('base', 'exist(''currentImage'', ''var'')')
        img = evalin('base', 'currentImage');
        app = guidata(gcf);
        app.OriginalImage = img;
        
        % 更新图像显示
        imshow(app.OriginalImage, 'Parent', app.OriginalAxes);
        title(app.OriginalAxes, '原始图像');
        
        % 应用当前方法
        applyEdgeDetection();
        
        guidata(gcf, app);
    end
end
