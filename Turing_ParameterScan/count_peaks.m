function nPeaks = count_peaks(Aend, x)
% 用 findpeaks 统计 A(x) 的峰数
% 可以根据需要调阈值，避免把小噪声算成峰

    % 粗略估一下典型高度，用中位数+一点 margin
    baseLevel = median(Aend);
    maxLevel  = max(Aend);
    % 设一个最低峰高阈值（介于 base 和 max 中间）
    minHeight = baseLevel + 0.2*(maxLevel - baseLevel);

    % 最小峰间距（物理上控制最小波长）
    minDist = (x(end)-x(1))/20;   % 这里允许最多 ~10 个峰，可根据需要改

    [pks, ~] = findpeaks(Aend, x, ...
                            'MinPeakHeight', minHeight, ...
                            'MinPeakDistance', minDist);

    nPeaks = numel(pks);
end
