%% check_peak5_candidates.m
clear; clc;

% load 参数候选
load('gm_peak5_candidates.mat','peak5_candidates','L');

nC = size(peak5_candidates,1);
true_peak5 = [];

for k = 1:nC
    aA  = peak5_candidates(k,1);
    aI  = peak5_candidates(k,2);
    rho = peak5_candidates(k,3);
    d   = peak5_candidates(k,4);

    fprintf('Candidate %2d / %2d: aA=%.3g, aI=%.3g, rho=%.3g, d=%.3g ... ', ...
            k, nC, aA, aI, rho, d);

    [x, Aend] = simulate_GM_1D(aA, aI, rho, d, L);

    nPeaks = count_peaks(Aend, x);

    fprintf('sim peaks = %d\n', nPeaks);

    if nPeaks == 5
        true_peak5 = [true_peak5; aA, aI, rho, d]; %#ok<AGROW>

        % 顺便画一张图
        figure;
        plot(x, Aend,'LineWidth',2);
        xlabel('x'); ylabel('A(x)');
        title(sprintf('5-peak pattern: a_A=%.3g, a_I=%.3g, \\rho=%.3g, d=%.3g', ...
              aA,aI,rho,d));
        grid on;
        drawnow;
    end
end

fprintf('最终真实模拟得到 5 个峰的参数共有 %d 组。\n', size(true_peak5,1));
% save('gm_true_peak5.mat','true_peak5','L');
