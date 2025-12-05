%% Part C: Spatial Sox9–Wnt pattern with A(x) & I(x) plotted together

clear; clc;

%% 空间离散
N  = 100;      % 空间点个数
L  = 1.0;      % 空间长度（0~1）
dx = L/(N-1);
x  = linspace(0, L, N);

%% 扩散系数（关键：DI >> DA）
DA = 1e-4;     % Sox9 diffusion
DI = 0.1;      % Wnt diffusion  (much larger)

%% toggle 参数（和 MatCont 一致）
alphaI = 4;    % Wnt最大生成率
m = 2;         % I 抑制 A 的 Hill 係数
n = 2;         % A 抑制 I 的 Hill 係数

% 这三个 alphaA 建议从 MatCont 分岔图里来，这里先给个示例
alphaA_WT     = 4.5;   % WT: 在双稳态区中间
alphaA_mut    = 6.0;   % Mutant: 右LP右侧，高A单稳态
alphaA_rescue = 5.0;   % Rescue: 比WT大一点，但仍在LP之间

alphaA_list   = [alphaA_WT, alphaA_mut, alphaA_rescue];
label_list    = {'WT','Mutant','Rescued'};

%% 初始条件：均匀+噪声
A0 = 0.5 + 0.1*randn(N,1);   % 噪声稍微大一点，容易长pattern
I0 = 0.5 + 0.1*randn(N,1);
Y0 = [A0; I0];

tspan = [0 1000];            % 时间要拉长一点，让图案稳定

%% 颜色设置
cA = [0.85 0.1 0.1];         % Sox9 颜色（红）
cI = [0.1 0.1 0.85];         % Wnt 颜色（蓝）

figure; clf;

for k = 1:3
    alphaA = alphaA_list(k);

    % 积分反应–扩散系统
    [t, Y] = ode15s(@(t,y) spatial_toggle_rhs(t,y,N,dx,DA,DI,alphaA,alphaI,m,n), ...
                    tspan, Y0);

    % 取最后一个时间点的空间分布
    A_final = Y(end, 1:N);
    I_final = Y(end, N+1:end);

    % 画图：同一 subplot 中画 Sox9 和 Wnt
    subplot(3,1,k);
    plot(x, A_final, '-', 'LineWidth', 2, 'Color', cA); hold on;
    plot(x, I_final, '--', 'LineWidth', 2, 'Color', cI);
    ylabel('Activity');
    title(label_list{k});
    if k == 1
        legend({'Sox9 (A)', 'Wnt (I)'}, 'Location','best');
    end
    if k == 3
        xlabel('Position along limb (x)');
    end
    box on;
end

sgtitle('Spatial Sox9 and Wnt patterns under WT, mutant, and rescued conditions');
