function edge_detection_gui()
% EDGE_DETECTION_GUI 边缘检测算法选择和参数调整界面
% 该GUI允许用户:
% 1. 加载图像
% 2. 选择边缘检测算法
% 3. 调整算法参数
% 4. 预览效果并应用到主程序

% 创建主界面
fig = figure('Name', '边缘检测算法选择', ...
    'Position', [200, 200, 1000, 600], ...
    'NumberTitle', 'off', ...
    'MenuBar', 'none', ...
    'ToolBar', 'none', ...
    'Resize', 'on', ...
    'CloseRequestFcn', @closeGUI);

% 全局变量
app.Fig = fig;
app.OriginalImage = [];
app.PreprocessedImage = [];
app.EdgeImage = [];
app.SelectedMethod = 'Canny';
app.CannyThreshLow = 0.1;
app.CannyThreshHigh = 0.2;
app.CannyStdDev = 1.0;
app.SobelThresh = 0.25;
app.PrewittThresh = 0.25;
app.RobertsThresh = 0.25;
app.LogThresh = 0.01;
app.LogStdDev = 2.0;
app.IsPreprocessed = true;
app.AppliedToMain = false;

% 创建面板
controlPanel = uipanel(fig, 'Title', '控制面板', ...
    'Position', [0.01, 0.01, 0.29, 0.98]);

imagePanel = uipanel(fig, 'Title', '图像预览', ...
    'Position', [0.31, 0.01, 0.68, 0.98]);

% 控制面板中的组件
% 图像加载按钮
uicontrol(controlPanel, 'Style', 'pushbutton', ...
    'String', '加载图像', ...
    'Position', [20, 540, 100, 30], ...
    'Callback', @loadImage);

% 退出按钮 - 添加退出按钮使用户可以关闭窗口
uicontrol(controlPanel, 'Style', 'pushbutton', ...
    'String', '退出', ...
    'Position', [130, 540, 100, 30], ...
    'BackgroundColor', [0.8 0.3 0.3], ...
    'Callback', @(~,~) delete(fig));

% 算法选择下拉框
uicontrol(controlPanel, 'Style', 'text', ...
    'String', '边缘检测算法:', ...
    'Position', [20, 500, 100, 20], ...
    'HorizontalAlignment', 'left');

app.MethodDropdown = uicontrol(controlPanel, 'Style', 'popupmenu', ...
    'String', {'Canny', 'Sobel', 'Prewitt', 'Roberts', 'LoG'}, ...
    'Position', [20, 470, 250, 25], ...
    'Callback', @methodChanged);

% 参数调整面板 - 初始显示Canny参数
app.ParamPanel = uipanel(controlPanel, 'Title', 'Canny 参数', ...
    'Position', [0.05, 0.3, 0.9, 0.28]);

% 预处理选择
app.PreprocessCheckbox = uicontrol(controlPanel, 'Style', 'checkbox', ...
    'String', '应用图像预处理', ...
    'Position', [20, 170, 150, 20], ...
    'Value', 1, ...
    'Callback', @togglePreprocess);

% 应用按钮
uicontrol(controlPanel, 'Style', 'pushbutton', ...
    'String', '应用到主程序', ...
    'Position', [75, 50, 120, 40], ...
    'BackgroundColor', [0.4 0.8 0.4], ...
    'Callback', @applyToMain);

% 比较所有算法按钮
uicontrol(controlPanel, 'Style', 'pushbutton', ...
    'String', '比较所有算法', ...
    'Position', [75, 100, 120, 30], ...
    'Callback', @compareAllMethods);

% 图像显示区域 - 创建更明确的全局对象
app.OriginalAxes = axes('Parent', imagePanel, 'Position', [0.05, 0.1, 0.42, 0.8]);
title(app.OriginalAxes, '原始图像');
axis(app.OriginalAxes, 'off');

app.EdgeAxes = axes('Parent', imagePanel, 'Position', [0.55, 0.1, 0.42, 0.8]);
title(app.EdgeAxes, '边缘检测效果');
axis(app.EdgeAxes, 'off');

% 初始化Canny参数控件
createCannyControls();

% 设置全局应用数据
guidata(fig, app);

% 回调函数
function loadImage(~, ~)
    app = guidata(fig);
    
    % 清空当前错误信息
    disp('尝试加载图像...');
    
    % 强制创建/更新轴对象 - 确保存在
    imagePanel = findobj(fig, 'Type', 'uipanel', 'Title', '图像预览');
    if ~isempty(imagePanel)
        if ~isfield(app, 'OriginalAxes') || isempty(app.OriginalAxes) || ~ishandle(app.OriginalAxes)
            % 创建新的轴对象
            app.OriginalAxes = axes('Parent', imagePanel, 'Position', [0.05, 0.1, 0.42, 0.8]);
            title(app.OriginalAxes, '原始图像');
        end
        
        if ~isfield(app, 'EdgeAxes') || isempty(app.EdgeAxes) || ~ishandle(app.EdgeAxes)
            app.EdgeAxes = axes('Parent', imagePanel, 'Position', [0.55, 0.1, 0.42, 0.8]);
            title(app.EdgeAxes, '边缘检测效果');
        end
    else
        errordlg('无法找到图像显示面板', '错误');
        return;
    end
    
    % 先尝试从工作区获取图像
    try
        if evalin('base', 'exist(''currentImage'', ''var'')')
            disp('从工作区加载图像');
            img = evalin('base', 'currentImage');
            if ~isempty(img)
                % 处理图像并显示
                processingAndDisplayImage(img);
                return;
            end
        end
    catch ME
        disp(['从工作区加载图像失败: ', ME.message]);
    end
    
    % 如果工作区没有图像，则打开文件选择对话框
    [filename, pathname] = uigetfile({'*.jpg;*.png;*.bmp;*.tif', '图像文件 (*.jpg, *.png, *.bmp, *.tif)'}, '选择一个图像文件');
    if isequal(filename, 0) || isequal(pathname, 0)
        return;
    end
    
    try
        img = imread(fullfile(pathname, filename));
        % 处理图像并显示
        processingAndDisplayImage(img);
    catch ME
        errordlg(['加载图像出错: ' ME.message], '错误');
    end
end

% 处理图像并显示 - 新增辅助函数
function processingAndDisplayImage(img)
    app = guidata(fig);
    
    % 确保图像是灰度图
    if size(img, 3) == 3
        img = rgb2gray(img);
    end
    
    app.OriginalImage = img;
    app.AppliedToMain = false;
    
    % 清空并重置轴对象
    if isfield(app, 'OriginalAxes') && ishandle(app.OriginalAxes)
        axes(app.OriginalAxes);
        cla;
        imshow(app.OriginalImage);
        title('原始图像');
    else
        disp('无法显示原始图像 - 轴对象不存在');
        return;
    end
    
    % 应用当前方法
    applyEdgeDetection();
    
    guidata(fig, app);
end

function methodChanged(~, ~)
    app = guidata(fig);
    methods = {'Canny', 'Sobel', 'Prewitt', 'Roberts', 'LoG'};
    app.SelectedMethod = methods{app.MethodDropdown.Value};
    app.AppliedToMain = false;
    
    % 删除旧参数控件
    delete(findall(app.ParamPanel, 'Type', 'uicontrol'));
    
    % 创建新参数控件
    switch app.SelectedMethod
        case 'Canny'
            app.ParamPanel.Title = 'Canny 参数';
            createCannyControls();
        case 'Sobel'
            app.ParamPanel.Title = 'Sobel 参数';
            createSobelControls();
        case 'Prewitt'
            app.ParamPanel.Title = 'Prewitt 参数';
            createPrewittControls();
        case 'Roberts'
            app.ParamPanel.Title = 'Roberts 参数';
            createRobertsControls();
        case 'LoG'
            app.ParamPanel.Title = 'LoG 参数';
            createLoGControls();
    end
    
    % 应用当前方法
    applyEdgeDetection();
    
    guidata(fig, app);
end

function createCannyControls()
    app = guidata(fig);
    
    % 低阈值
    uicontrol(app.ParamPanel, 'Style', 'text', ...
        'String', '低阈值:', ...
        'Position', [10, 120, 80, 20], ...
        'HorizontalAlignment', 'left');
    
    app.CannyLowSlider = uicontrol(app.ParamPanel, 'Style', 'slider', ...
        'Min', 0, 'Max', 1, 'Value', app.CannyThreshLow, ...
        'Position', [100, 120, 100, 20], ...
        'Callback', @(hObject,~) updateCannyLowValue(hObject.Value));
    
    app.CannyLowText = uicontrol(app.ParamPanel, 'Style', 'text', ...
        'String', num2str(app.CannyThreshLow, '%.2f'), ...
        'Position', [210, 120, 40, 20]);
    
    % 高阈值
    uicontrol(app.ParamPanel, 'Style', 'text', ...
        'String', '高阈值:', ...
        'Position', [10, 90, 80, 20], ...
        'HorizontalAlignment', 'left');
    
    app.CannyHighSlider = uicontrol(app.ParamPanel, 'Style', 'slider', ...
        'Min', 0, 'Max', 1, 'Value', app.CannyThreshHigh, ...
        'Position', [100, 90, 100, 20], ...
        'Callback', @(hObject,~) updateCannyHighValue(hObject.Value));
    
    app.CannyHighText = uicontrol(app.ParamPanel, 'Style', 'text', ...
        'String', num2str(app.CannyThreshHigh, '%.2f'), ...
        'Position', [210, 90, 40, 20]);
    
    % Sigma (标准差)
    uicontrol(app.ParamPanel, 'Style', 'text', ...
        'String', 'Sigma:', ...
        'Position', [10, 60, 80, 20], ...
        'HorizontalAlignment', 'left');
    
    app.CannySigmaSlider = uicontrol(app.ParamPanel, 'Style', 'slider', ...
        'Min', 0.1, 'Max', 5, 'Value', app.CannyStdDev, ...
        'Position', [100, 60, 100, 20], ...
        'Callback', @(hObject,~) updateCannySigmaValue(hObject.Value));
    
    app.CannySigmaText = uicontrol(app.ParamPanel, 'Style', 'text', ...
        'String', num2str(app.CannyStdDev, '%.2f'), ...
        'Position', [210, 60, 40, 20]);
    
    guidata(fig, app);
end

function createSobelControls()
    app = guidata(fig);
    
    % 阈值
    uicontrol(app.ParamPanel, 'Style', 'text', ...
        'String', '阈值:', ...
        'Position', [10, 90, 80, 20], ...
        'HorizontalAlignment', 'left');
    
    app.SobelSlider = uicontrol(app.ParamPanel, 'Style', 'slider', ...
        'Min', 0, 'Max', 1, 'Value', app.SobelThresh, ...
        'Position', [100, 90, 100, 20], ...
        'Callback', @(hObject,~) updateSobelValue(hObject.Value));
    
    app.SobelText = uicontrol(app.ParamPanel, 'Style', 'text', ...
        'String', num2str(app.SobelThresh, '%.2f'), ...
        'Position', [210, 90, 40, 20]);
    
    guidata(fig, app);
end

function createPrewittControls()
    app = guidata(fig);
    
    % 阈值
    uicontrol(app.ParamPanel, 'Style', 'text', ...
        'String', '阈值:', ...
        'Position', [10, 90, 80, 20], ...
        'HorizontalAlignment', 'left');
    
    app.PrewittSlider = uicontrol(app.ParamPanel, 'Style', 'slider', ...
        'Min', 0, 'Max', 1, 'Value', app.PrewittThresh, ...
        'Position', [100, 90, 100, 20], ...
        'Callback', @(hObject,~) updatePrewittValue(hObject.Value));
    
    app.PrewittText = uicontrol(app.ParamPanel, 'Style', 'text', ...
        'String', num2str(app.PrewittThresh, '%.2f'), ...
        'Position', [210, 90, 40, 20]);
    
    guidata(fig, app);
end

function createRobertsControls()
    app = guidata(fig);
    
    % 阈值
    uicontrol(app.ParamPanel, 'Style', 'text', ...
        'String', '阈值:', ...
        'Position', [10, 90, 80, 20], ...
        'HorizontalAlignment', 'left');
    
    app.RobertsSlider = uicontrol(app.ParamPanel, 'Style', 'slider', ...
        'Min', 0, 'Max', 1, 'Value', app.RobertsThresh, ...
        'Position', [100, 90, 100, 20], ...
        'Callback', @(hObject,~) updateRobertsValue(hObject.Value));
    
    app.RobertsText = uicontrol(app.ParamPanel, 'Style', 'text', ...
        'String', num2str(app.RobertsThresh, '%.2f'), ...
        'Position', [210, 90, 40, 20]);
    
    guidata(fig, app);
end

function createLoGControls()
    app = guidata(fig);
    
    % 阈值
    uicontrol(app.ParamPanel, 'Style', 'text', ...
        'String', '阈值:', ...
        'Position', [10, 90, 80, 20], ...
        'HorizontalAlignment', 'left');
    
    app.LogSlider = uicontrol(app.ParamPanel, 'Style', 'slider', ...
        'Min', 0, 'Max', 0.1, 'Value', app.LogThresh, ...
        'Position', [100, 90, 100, 20], ...
        'Callback', @(hObject,~) updateLoGThreshValue(hObject.Value));
    
    app.LogText = uicontrol(app.ParamPanel, 'Style', 'text', ...
        'String', num2str(app.LogThresh, '%.4f'), ...
        'Position', [210, 90, 50, 20]);
    
    % Sigma (标准差)
    uicontrol(app.ParamPanel, 'Style', 'text', ...
        'String', 'Sigma:', ...
        'Position', [10, 60, 80, 20], ...
        'HorizontalAlignment', 'left');
    
    app.LogSigmaSlider = uicontrol(app.ParamPanel, 'Style', 'slider', ...
        'Min', 0.1, 'Max', 5, 'Value', app.LogStdDev, ...
        'Position', [100, 60, 100, 20], ...
        'Callback', @(hObject,~) updateLoGSigmaValue(hObject.Value));
    
    app.LogSigmaText = uicontrol(app.ParamPanel, 'Style', 'text', ...
        'String', num2str(app.LogStdDev, '%.2f'), ...
        'Position', [210, 60, 40, 20]);
    
    guidata(fig, app);
end

function updateCannyLowValue(value)
    app = guidata(fig);
    app.CannyThreshLow = value;
    app.CannyLowText.String = num2str(value, '%.2f');
    app.AppliedToMain = false;
    applyEdgeDetection();
    guidata(fig, app);
end

function updateCannyHighValue(value)
    app = guidata(fig);
    app.CannyThreshHigh = value;
    app.CannyHighText.String = num2str(value, '%.2f');
    app.AppliedToMain = false;
    applyEdgeDetection();
    guidata(fig, app);
end

function updateCannySigmaValue(value)
    app = guidata(fig);
    app.CannyStdDev = value;
    app.CannySigmaText.String = num2str(value, '%.2f');
    app.AppliedToMain = false;
    applyEdgeDetection();
    guidata(fig, app);
end

function updateSobelValue(value)
    app = guidata(fig);
    app.SobelThresh = value;
    app.SobelText.String = num2str(value, '%.2f');
    app.AppliedToMain = false;
    applyEdgeDetection();
    guidata(fig, app);
end

function updatePrewittValue(value)
    app = guidata(fig);
    app.PrewittThresh = value;
    app.PrewittText.String = num2str(value, '%.2f');
    app.AppliedToMain = false;
    applyEdgeDetection();
    guidata(fig, app);
end

function updateRobertsValue(value)
    app = guidata(fig);
    app.RobertsThresh = value;
    app.RobertsText.String = num2str(value, '%.2f');
    app.AppliedToMain = false;
    applyEdgeDetection();
    guidata(fig, app);
end

function updateLoGThreshValue(value)
    app = guidata(fig);
    app.LogThresh = value;
    app.LogText.String = num2str(value, '%.4f');
    app.AppliedToMain = false;
    applyEdgeDetection();
    guidata(fig, app);
end

function updateLoGSigmaValue(value)
    app = guidata(fig);
    app.LogStdDev = value;
    app.LogSigmaText.String = num2str(value, '%.2f');
    app.AppliedToMain = false;
    applyEdgeDetection();
    guidata(fig, app);
end

function togglePreprocess(hObject, ~)
    app = guidata(fig);
    app.IsPreprocessed = hObject.Value;
    app.AppliedToMain = false;
    applyEdgeDetection();
    guidata(fig, app);
end

function applyEdgeDetection()
    app = guidata(fig);
    
    if isempty(app.OriginalImage)
        return;
    end
    
    % 检查轴对象是否存在
    if ~isfield(app, 'EdgeAxes') || ~ishandle(app.EdgeAxes)
        disp('边缘图像显示区域不存在，无法显示结果');
        return;
    end
    
    % 预处理图像
    if app.IsPreprocessed
        % 应用更强的预处理来减少噪声
        preprocessed_img = imgaussfilt(app.OriginalImage, 1.5);
        % 提高对比度
        preprocessed_img = imadjust(preprocessed_img);
        % 中值滤波消除椒盐噪声
        preprocessed_img = medfilt2(preprocessed_img, [3 3]);
    else
        preprocessed_img = app.OriginalImage;
    end
    
    app.PreprocessedImage = preprocessed_img;
    
    try
        % 应用边缘检测
        switch app.SelectedMethod
            case 'Canny'
                edges = edge(preprocessed_img, 'Canny', [app.CannyThreshLow app.CannyThreshHigh], app.CannyStdDev);
            case 'Sobel'
                edges = edge(preprocessed_img, 'Sobel', app.SobelThresh);
            case 'Prewitt'
                edges = edge(preprocessed_img, 'Prewitt', app.PrewittThresh);
            case 'Roberts'
                edges = edge(preprocessed_img, 'Roberts', app.RobertsThresh);
            case 'LoG'
                edges = edge(preprocessed_img, 'log', app.LogThresh, app.LogStdDev);
        end
        
        app.EdgeImage = edges;
        
        % 显示结果 - 重置轴并显示
        axes(app.EdgeAxes);
        cla(app.EdgeAxes); % 清除当前轴内容
        imshow(app.EdgeImage, 'Parent', app.EdgeAxes);
        title(app.EdgeAxes, [app.SelectedMethod, ' 边缘检测结果']);
        
        guidata(fig, app);
    catch ME
        disp(['边缘检测错误: ', ME.message]);
        errordlg(['边缘检测错误: ', ME.message], '错误');
    end
end

function compareAllMethods(~, ~)
    app = guidata(fig);
    
    if isempty(app.OriginalImage)
        warndlg('请先加载图像', '警告');
        return;
    end
    
    % 调用比较函数
    compare_edge_detectors(app.OriginalImage);
end

function applyToMain(~, ~)
    app = guidata(fig);
    
    if isempty(app.EdgeImage)
        warndlg('请先加载图像并检测边缘', '警告');
        return;
    end
    
    % 创建结构体存储所选方法和参数
    edgeInfo = struct();
    edgeInfo.method = app.SelectedMethod;
    edgeInfo.image = app.EdgeImage;
    edgeInfo.isPreprocessed = app.IsPreprocessed;
    
    switch app.SelectedMethod
        case 'Canny'
            edgeInfo.params = struct('low', app.CannyThreshLow, 'high', app.CannyThreshHigh, 'sigma', app.CannyStdDev);
        case 'Sobel'
            edgeInfo.params = struct('threshold', app.SobelThresh);
        case 'Prewitt'
            edgeInfo.params = struct('threshold', app.PrewittThresh);
        case 'Roberts'
            edgeInfo.params = struct('threshold', app.RobertsThresh);
        case 'LoG'
            edgeInfo.params = struct('threshold', app.LogThresh, 'sigma', app.LogStdDev);
    end
    
    % 保存边缘检测结果到工作区，使主程序可以访问
    assignin('base', 'edgeDetectionResult', edgeInfo);
    
    app.AppliedToMain = true;
    guidata(fig, app);
    
    % 提示用户结果已应用
    msgbox(['已将 ', app.SelectedMethod, ' 边缘检测结果应用到主程序。请在主程序中通过 edgeDetectionResult 变量访问结果。'], '成功');
end

function closeGUI(~, ~)
    app = guidata(fig);
    
    % 直接删除图形窗口，不再弹出确认对话框，简化退出流程
    try
        % 如果有结果但未应用，自动应用
        if ~app.AppliedToMain && ~isempty(app.EdgeImage)
            applyToMain([], []);
        end
    catch
        % 忽略错误，确保窗口可以关闭
    end
    
    % 直接删除图形窗口
    try
        delete(fig);
    catch
        % 如果常规删除失败，尝试强制关闭
        if ishandle(fig)
            close(fig, 'force');
        end
    end
end

end
