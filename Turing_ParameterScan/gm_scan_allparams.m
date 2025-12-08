%% gm_scan_allparams.m
% 扫描无量纲 GM 模型参数，寻找图灵斑纹参数集

clear; clc;

%----- 仿真设置 -----
L     = 10.0;   % 无量纲空间长度，用于确定可用的空间模态 k_n = n*pi/L
nMaxK = 40;    % 扫描的最大 mode index (n=1..nMaxK)

Nsample = 30000;   % 随机采样的参数个数，可酌情调大/调小

%----- 参数采样范围（无量纲） -----
% 采用 log-uniform 采样，让多个数量级都被覆盖
aA_min = 0.1;   aA_max = 20;     % aA = alpha_A / mu_A
aI_min = 0.1;   aI_max = 20;     % aI = alpha_I / mu_A
rho_min = 1.01; rho_max = 20;    % rho = mu_I / mu_A (>1 才可能图灵)
d_min  = 1;   d_max  = 100;    % d = D_I / D_A (扩散比)

params  = zeros(Nsample,4);   % 每行存 [aA, aI, rho, d]
isTuring = false(Nsample,1);  % 是否满足图灵不稳定

fprintf('Start sampling %d parameter sets...\n', Nsample);

for n = 1:Nsample
    % log-uniform 抽样
    aA  = 10^( log10(aA_min) + (log10(aA_max)-log10(aA_min))*rand );
    aI  = 10^( log10(aI_min) + (log10(aI_max)-log10(aI_min))*rand );
    rho = 10^( log10(rho_min)+ (log10(rho_max)-log10(rho_min))*rand );
    d   = 10^( log10(d_min) + (log10(d_max)-log10(d_min))*rand );

    params(n,:) = [aA, aI, rho, d];

    % 检查该组参数是否产生图灵不稳定
    isTuring(n) = isTuringGM(aA, aI, rho, d, L, nMaxK);
end

fprintf('Done. Turing points: %d / %d\n', sum(isTuring), Nsample);

%% 2D 投影：rho - d 平面上的图灵区域
figure;
hold on;
% 先画非图灵点（灰色）
scatter(params(~isTuring,3), params(~isTuring,4), 10, [0.8 0.8 0.8], 'filled');
% 再画图灵点（红色）
scatter(params(isTuring,3), params(isTuring,4), 10, 'r', 'filled');

set(gca,'XScale','log','YScale','log');
xlabel('\rho = \mu_I / \mu_A');
ylabel('d = D_I / D_A');
title('Gierer–Meinhardt Turing region (projection onto \rho–d plane)');
legend('no pattern','Turing','Location','best');
grid on;
hold off;

% 为了方便，先取索引
idxT  = isTuring;
idxNT = ~isTuring;

%% 1) aA - aI 平面
figure;
hold on;
scatter(params(idxNT,1), params(idxNT,2), 10, [0.8 0.8 0.8], 'filled');
scatter(params(idxT ,1), params(idxT ,2), 10, 'r', 'filled');
set(gca,'XScale','log','YScale','log');
xlabel('\alpha_A / \mu_A');
ylabel('\alpha_I / \mu_A');
title('Projection onto a_A - a_I plane');
legend('no pattern','Turing','Location','best');
grid on;
hold off;

%% 2) aA - d 平面
figure;
hold on;
scatter(params(idxNT,1), params(idxNT,4), 10, [0.8 0.8 0.8], 'filled');
scatter(params(idxT ,1), params(idxT ,4), 10, 'r', 'filled');
set(gca,'XScale','log','YScale','log');
xlabel('\alpha_A / \mu_A');
ylabel('d = D_I / D_A');
title('Projection onto a_A - d plane');
legend('no pattern','Turing','Location','best');
grid on;
hold off;

%% 3) aI - d 平面
figure;
hold on;
scatter(params(idxNT,2), params(idxNT,4), 10, [0.8 0.8 0.8], 'filled');
scatter(params(idxT ,2), params(idxT ,4), 10, 'r', 'filled');
set(gca,'XScale','log','YScale','log');
xlabel('\alpha_I / \mu_A');
ylabel('d = D_I / D_A');
title('Projection onto a_I - d plane');
legend('no pattern','Turing','Location','best');
grid on;
hold off;



figure;
scatter3(params(idxNT,1), params(idxNT,2), params(idxNT,4), ...
         10, [0.8 0.8 0.8], 'filled');  % 先画灰色
hold on;
scatter3(params(idxT ,1), params(idxT ,2), params(idxT ,4), ...
         10, 'r', 'filled');            % 再画红色
set(gca,'XScale','log','YScale','log','ZScale','log');
xlabel('\alpha_A / \mu_A');
ylabel('\alpha_I / \mu_A');
zlabel('d = D_I / D_A');
title('3D view: (a_A, a_I, d)');
grid on;
view(45,30);   % 调一下视角
legend('no pattern','Turing','Location','best');
hold off;

save('gm_param_scan.mat', 'params', 'isTuring');
