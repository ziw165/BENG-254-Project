%%find_peak5_candidates.m
clear; clc;

% 如果你之前把结果存成 mat 文件了：
load('gm_param_scan.mat','params','isTuring');

if ~exist('params','var') || ~exist('isTuring','var')
    error('请先在 workspace 中准备 params 和 isTuring');
end

L = 40;              % 你 PDE 模拟时的空间长度，要跟 simulate_GM_1D 里面一致
targetPeaks = 5;

idxT = find(isTuring);
nT   = numel(idxT);

peak5_candidates = [];

fprintf('从 %d 个 Turing 参数中预测峰数...\n', nT);

for k = 1:nT
    p = params(idxT(k),:);
    aA  = p(1);
    aI  = p(2);
    rho = p(3);
    d   = p(4);

    nPred = predict_peaks_GM(aA,aI,rho,d,L);

    if nPred == targetPeaks
        peak5_candidates = [peak5_candidates; p]; %#ok<AGROW>
    end
end

fprintf('预测峰数 = %d 的候选参数有 %d 组。\n', ...
        targetPeaks, size(peak5_candidates,1));

% 可选：保存下来
% save('gm_peak5_candidates.mat','peak5_candidates','L');
