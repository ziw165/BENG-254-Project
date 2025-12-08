%% gm_2D_pattern.m
% 2D Turing pattern of the dimensionless Gierer–Meinhardt model
% dA/dt = aA*A^2/I - A + DA * Lap(A)
% dI/dt = aI*A^2   - rho*I + DI * Lap(I)

clear; clc; close all;

%% ---------- 1. Parameters ----------

% reaction parameters (choose one parameter set in the Turing region)
aA  = 0.178;    % alpha_A / mu_A
aI  = 8.86;     % alpha_I / mu_A
rho = 3.0;      % mu_I  / mu_A
d   = 40;       % D_I / D_A  (moderate value for 2D demo; can increase with smaller dt)

DA = 1;         % activator diffusion (dimensionless)
DI = d;         % inhibitor diffusion

% domain and grid
Lx = 60;        % domain size in x
Ly = 40;        % domain size in y
Nx = 128;       % number of grid points in x
Ny = 128;       % number of grid points in y

dx = Lx / Nx;
dy = Ly / Ny;

x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);

% time stepping
dt    = 1e-4;   % time step (may need smaller if DI is very large)
Tend  = 15;     % final time
nStep = round(Tend / dt);
plotEvery = 500;  % how often to update figure

fprintf('Running 2D GM simulation: %d x %d grid, %d steps...\n', Nx, Ny, nStep);

%% ---------- 2. Initial condition: steady state + noise ----------
% ----- 初始条件：稳态 + 一维余弦扰动 → 自然长成条纹 -----
Astar = aA * rho / aI;
Istar = aA^2 * rho / aI;

[X,Y] = meshgrid(x,y);

eps0 = 0.05;   % 条纹扰动幅度，5% 左右
k = 2*pi / (Lx/5);   % 控制条纹条数，这里大约 5 条

% 载入 1D 模式
load('A1D.mat');   % 假设你把 A1D 保存下来了
load('I1D.mat');   % 同理保存 I1D

% 或者如果你直接在同一工作区运行，可以跳过 load，直接用 A1D

% 将 1D pattern 延伸到 2D
A = repmat(A1D, Ny, 1);   % 每一行都是 A1D
I = repmat(I1D, Ny, 1);   % 每一行都是 I1D

% 加一点非常小的噪声（可选）
A = A .* (1 + 0.001*randn(size(A)));
I = I .* (1 + 0.001*randn(size(I)));



%% ---------- 3. Time integration (explicit Euler, periodic BC) ----------

figure('Position',[100 100 800 700]);

for n = 1:nStep

    % periodic Laplacian via circshift
    LapA = ( circshift(A,[0  1]) + circshift(A,[0 -1]) + ...
             circshift(A,[ 1 0]) + circshift(A,[-1 0]) - 4*A ) / dx^2;

    LapI = ( circshift(I,[0  1]) + circshift(I,[0 -1]) + ...
             circshift(I,[ 1 0]) + circshift(I,[-1 0]) - 4*I ) / dx^2;

    % reaction terms
    % (avoid division by zero)
    I_safe = max(I, 1e-8);

    R_A = aA * A.^2 ./ I_safe - A;
    R_I = aI * A.^2           - rho*I;

    % explicit Euler update
    A = A + dt*( R_A + DA*LapA );
    I = I + dt*( R_I + DI*LapI );

    % keep concentrations non-negative
    A = max(A, 0);
    I = max(I, 0);

    % visualization
    if mod(n, plotEvery) == 0 || n == 1
        imagesc(x, y, A);
        axis image; axis xy;
        colorbar;
        clim([0, max(Astar*2, max(A(:)))]);
        title(sprintf('Activator A(x,y), t = %.2f', n*dt));
        xlabel('x'); ylabel('y');
        drawnow;
    end
end

fprintf('Done.\n');
