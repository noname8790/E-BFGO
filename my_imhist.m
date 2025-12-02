function [counts, x] = my_imhist(I, n)
% MY_IMHIST 自定义图像直方图计算函数
%   参数:
%       I - 输入图像(灰度或RGB)
%       n - 直方图bin数量(可选，默认256)
%   返回值:
%       counts - 直方图统计值
%       binLocations - bin中心位置

if nargin < 2
    n = 256; % 默认256个bin
end

if ismatrix(I) % 灰度图像
    if isinteger(I)
        range = [intmin(class(I)) intmax(class(I))];
    else
        range = [0 1]; % 归一化图像
    end
    
    [counts, edges] = histcounts(I(:), n, 'BinLimits', range);
    % x = (edges(1:end-1) + edges(2:end)) / 2;
    
elseif ndims(I) == 3 && size(I,3) == 3 % RGB图像
    counts = zeros(n, 3);
    for ch = 1:3
        [counts(:,ch), edges] = histcounts(I(:,:,ch), n);
        % if ch == 1
        %     x = (edges(1:end-1) + edges(2:end)) / 2;
        % end
    end
else
    error('不支持的图像格式');
end
x = 0:255;
end