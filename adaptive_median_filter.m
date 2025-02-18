function filtered = adaptive_median_filter(img)
% 自适应中值滤波实现
[m, n] = size(img);
filtered = zeros(m, n);
min_window = 3;
max_window = 9;  % 增加最大窗口大小

for i = 1:m
    for j = 1:n
        window_size = min_window;
        while window_size <= max_window
            window = get_window(img, i, j, window_size);
            z_min = min(window(:));
            z_max = max(window(:));
            z_med = median(window(:));
            z_xy = img(i,j);
            
            if z_med > z_min && z_med < z_max
                if z_xy > z_min && z_xy < z_max
                    filtered(i,j) = z_xy;
                else
                    filtered(i,j) = z_med;
                end
                break;
            else
                window_size = window_size + 2;
                if window_size > max_window
                    filtered(i,j) = z_med;
                end
            end
        end
    end
end

% 增强对比度
filtered = imadjust(filtered);
filtered = uint8(filtered);
end

function window = get_window(img, i, j, window_size)
% 获取图像窗口
half = floor(window_size/2);
[m, n] = size(img);
row_min = max(1, i-half);
row_max = min(m, i+half);
col_min = max(1, j-half);
col_max = min(n, j+half);
window = img(row_min:row_max, col_min:col_max);
end 