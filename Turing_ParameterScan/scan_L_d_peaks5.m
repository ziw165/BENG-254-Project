%% scan_L_d_peak5.m
% 扫描 L (系统长度) 与 d = D_I/D_A (扩散比)
% 统计图灵模式的峰数，并画出 peaks = 5 的参数域

clear; clc;

%% -------- Step 1: 固定其它参数（典型 Turing 点）---------
% 建议你填入自己找到的那组产生 Turing 的参数：
aA  = 0.178;      % alpha_A / mu_A
aI  = 8.863;      % alpha_I / mu_A
rho = 3.0;        % mu_I / mu_A   （确保 >1）
% 你可以随时修改这些值，看不同参数下的 L-d 域变化

%% -------- Step 2: 扫描范围设置 ---------
L_list = linspace(20,150,80);      % 扫不同长度，例如 20~80
d_list = logspace(1, 3, 80);     % d = 0.1 ~ 1000 的 log 扫描

nL = length(L_list);
nd = length(d_list);

peakMap = nan(nL, nd);            % 保存峰数

fprintf('开始扫描 L × d 空间...\n');

%% -------- Step 3: 双循环扫描 ---------
for i = 1:nL
    L = L_list(i);
    for j = 1:nd
        d = d_list(j);

        %---- PDE 模拟 ----%
        [x, Aend] = simulate_GM_1D(aA, aI, rho, d, L);

        %---- 数峰 ----%
        nPeaks = count_peaks(Aend, x);

        peakMap(i,j) = nPeaks;

        fprintf('L = %.1f, d = %.3g → peaks = %d\n', L, d, nPeaks);
    end
end

fprintf('扫描结束。\n');

%% -------- Step 4: 绘制二维参数域图 ---------
figure('Position',[100 100 1600 1200]);
hold on;

% ---- 已经有的部分：画热图 + 红点 ----
imagesc(L_list, d_list, peakMap');    
set(gca,'YScale','log','YDir','normal');
colorbar;
xlabel('L (system length)');
ylabel('d = D_I / D_A');
title('Peak number map in (L, d) space');
hold on;

[row5, col5] = find(peakMap == 5);
L5 = L_list(row5);
d5 = d_list(col5);

scatter(L5, d5, 20, 'r', 'filled');
legend('peak count','peaks = 5');

% ========= 新增部分：拟合红点成曲线 =========

% 1) 在 log10(d) 空间里拟合二次多项式：log10(d) = a0 + a1*L + a2*L^2
logd5 = log10(d5);
p = polyfit(L5(:), logd5(:), 2);     % poly2 你也可以改成 1 或 3 看效果

% 2) 生成平滑的 L 取值
L_fit = linspace(min(L5), max(L5), 300);

% 3) 计算拟合曲线对应的 d
logd_fit = polyval(p, L_fit);
d_fit    = 10.^logd_fit;

% 4) 画在图上
plot(L_fit, d_fit, 'w-', 'LineWidth', 2);   % 白色粗线
% 也可以改成黑色: 'k-'

% （可选）在命令行看看拟合系数
disp('polyfit coefficients for log10(d) = a0 + a1*L + a2*L^2:');
disp(p);
