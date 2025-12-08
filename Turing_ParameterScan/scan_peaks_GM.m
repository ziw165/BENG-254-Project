%% scan_peaks_GM.m
% 在已知的 Turing 参数集中，筛选出 1D 空间模式 peak 数 = 5 的参数

clear; clc;

% 如果你是单独运行这个脚本，可以把原来保存的 params / isTuring 加载进来：
load('gm_param_scan.mat','params','isTuring');

% 这里假设 params 和 isTuring 已经在 workspace 中
if ~exist('params','var') || ~exist('isTuring','var')
    error('请先在工作区中生成 params 和 isTuring（4D 扫描结果）');
end

idxT = find(isTuring);          % 所有 Turing 参数的索引
nT   = numel(idxT);

% 为了不太慢，只挑一部分 Turing 点来模拟（可以改大/改小）
maxSim = 300;                   % 最多模拟多少个参数点
nSim   = min(nT, maxSim);

% 随机挑 nSim 个 Turing 点（也可以用 1:nSim 取前几个）
rng(1);                         % 固定随机种子，结果可复现
selIdx = idxT(randperm(nT, nSim));

targetPeaks = 5;               % 目标峰数

peakParams = [];               % 存 peak 数 = 5 的参数
peakInfo   = [];               % 额外信息（比如波长等，可以后扩展）

fprintf('将从 %d 个 Turing 参数中抽取 %d 个进行 1D PDE 模拟...\n', nT, nSim);

for k = 1:nSim
    p = params(selIdx(k),:);
    aA  = p(1);
    aI  = p(2);
    rho = p(3);
    d   = p(4);

    fprintf('Sim %3d / %3d: aA=%.3g, aI=%.3g, rho=%.3g, d=%.3g ... ', ...
            k, nSim, aA, aI, rho, d);

    %---- 用子函数模拟 1D GM 模型 ----%
    [x, Aend] = simulate_GM_1D(aA, aI, rho, d);

    %---- 统计峰数 ----%
    nPeaks = count_peaks(Aend, x);

    fprintf('peaks = %d\n', nPeaks);

    if nPeaks == targetPeaks
        peakParams = [peakParams; p];
        peakInfo   = [peakInfo; nPeaks]; %#ok<AGROW>

        % 可以顺便画一下这个模式（可选）
        figure;
        plot(x, Aend, 'LineWidth', 2);
        xlabel('x'); ylabel('A(x)');
        title(sprintf('A(x) with %d peaks (a_A=%.3g, a_I=%.3g, \\rho=%.3g, d=%.3g)', ...
              nPeaks, aA, aI, rho, d));
        grid on;
        drawnow;
    end
end

fprintf('\n筛选完成：共有 %d 组参数的 peak 数 = %d。\n', size(peakParams,1), targetPeaks);

% 想保存下来可以：
% save('gm_peak5_params.mat','peakParams','peakInfo');
