%% gm_1D_pattern_demo.m
% 无量纲 Gierer–Meinhardt 模型在 1D 上的图灵斑纹示例
% dA/dt = aA*A^2/I - A + d2A/dx2
% dI/dt = aI*A^2   - rho*I + d*d2I/dx2

clear; clc;

%% ===== 1. 参数设置（把这里改成你选中的那一点） =====
aA  = 0.178;    % alpha_A / mu_A  (X 轴)
aI  = 8.863;    % alpha_I / mu_A  (Y 轴)
rho = 3.0;      % mu_I / mu_A     （你需要填成对应的值，只要 >1 即可出图灵）
d    = 94;   % D_I / D_A       (Z 轴)

% 空间与时间
L  = 60;              % 1D 空间长度（无量纲，取大一点可以容纳多个峰）
Nx = 200;             % 空间网格数
x  = linspace(0,L,Nx);

Tend = 300;           % 终止时间（无量纲）
Nt   = 100;           % 采样的时间点数量（用于保存结果）
tspan = linspace(0, Tend, Nt);

% 初始扰动大小（相对于均一稳态）
noise_amp = 0.01;     % 1% 随机扰动

%% ===== 2. 计算均一稳态，用作初值中心 =====
Astar = aA * rho / aI;
Istar = aA^2 * rho / aI;

fprintf('Steady state: A* = %.4f, I* = %.4f\n', Astar, Istar);

%% ===== 3. 调用 pdepe 求解 1D PDE =====
% u = [A; I]

m = 0;  % 1D 笛卡尔坐标

sol = pdepe(m, ...
    @(x,t,u,dudx) gm_pdefun(x,t,u,dudx,aA,aI,rho,d), ...
    @(x)          gm_icfun(x,Astar,Istar,noise_amp), ...
    @(xl,ul,xr,ur,t) gm_bcfun(xl,ul,xr,ur,t), ...
    x, tspan);

% sol 的维度: [Nt, Nx, 2]
A = squeeze(sol(:,:,1));
I = squeeze(sol(:,:,2));

%% ===== 4. 画最终时刻的空间浓度曲线 =====
figure;
plot(x, A(end,:), 'LineWidth', 2); hold on;
plot(x, I(end,:), '--', 'LineWidth', 2);
xlabel('x (space)');
ylabel('Concentration');
legend('Activator A','Inhibitor I');
title(sprintf('GM pattern at t = %.1f (a_A=%.3g, a_I=%.3g, \\rho=%.3g, d=%.3g)', ...
      tspan(end), aA, aI, rho, d));
grid on;

%% ===== 5. 可选：画 A(x,t) 的时空图像 =====
figure;
imagesc(x, tspan, A); axis xy;
xlabel('x'); ylabel('t');
title('Activator A(x,t)');
colorbar;

%% ====== 局部函数定义 ======
function [c,f,s] = gm_pdefun(~,~,u,dudx,aA,aI,rho,d)
    % u(1) = A, u(2) = I
    A = u(1);
    I = u(2);
    Ax = dudx(1);
    Ix = dudx(2);

    % c(x,t,u, dudx) * ut = d/dx( f ) + s
    c = [1; 1];

    % 扩散通量 f
    f = [Ax; d*Ix];   % D_A=1, D_I=d

    % 反应项 s
    s1 = aA*A^2 / I - A;
    s2 = aI*A^2     - rho*I;
    s  = [s1; s2];
end

function u0 = gm_icfun(~,Astar,Istar,noise_amp)
    % 初始条件：稳态附近加一点空间噪声
    A0 = Astar * (1 + noise_amp * (2*rand-1));
    I0 = Istar * (1 + noise_amp * (2*rand-1));
    u0 = [A0; I0];
end

function [pl,ql,pr,qr] = gm_bcfun(~,~,~,~,~)
    % Neumann 边界：A_x = 0, I_x = 0
    % 在 pdepe 里，用 flux = 0 表示：p = 0, q = 1
    pl = [0;0];
    ql = [1;1];
    pr = [0;0];
    qr = [1;1];
end

A1D = A(end, :);
I1D = I(end, :);
save('A1D.mat', 'A1D');
save('I1D.mat', 'I1D');
