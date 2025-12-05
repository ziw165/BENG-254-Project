%% Part C: Turing-like digit pattern in 1D using activator-inhibitor RD

 close all;

%% Spatial domain
N = 600;          % number of spatial points
L = 3.9;          % limb length
dx = L/(N-1);
x = linspace(0, L, N);

%% Time
tspan = [0 500];

%% Parameters
alphaA = 4;
alphaI = 4;
muA    = 1;
muI    = 1;

DA = 1e-3;        % slow activator diffusion
DI = 1e-2;        % fast inhibitor diffusion

%% Initial conditions (small random noise)
A0 = 1e-3 * ones(N,1) + 1e-4 * cos(2*pi*x/L)';
I0 = 1e-3 * ones(N,1) + 1e-4 * cos(2*pi*x/L)';

Y0 = [A0; I0];

%% Solve PDE via method of lines
[t,Y] = ode15s(@(t,y) rd_rhs(t,y,N,dx,DA,DI,alphaA,alphaI,muA,muI), tspan, Y0);

%% Final steady state pattern
A_final = Y(end,1:N);
I_final = Y(end,N+1:end);

%% Plot spatial pattern
figure;
plot(x, A_final, 'r', 'LineWidth',2); hold on;
plot(x, I_final, 'b--', 'LineWidth',2);
xlabel('Position along limb (x)');
ylabel('Activity');
legend('Sox9-like activator A','Wnt-like inhibitor I');
title('Turing-like digit pattern (Part C)');
grid on;
%% Count number of Sox9 peaks (digit condensations)
[~, locs] = findpeaks(A_final, ...
    'MinPeakProminence', 0.05, ...        % 峰至少要比周围高这么多
    'MinPeakDistance',  round(0.1/dx));   % 相邻峰之间至少隔 0.1 空间单位

nPeaks = numel(locs);
fprintf('Number of Sox9 peaks (digits) = %d\n', nPeaks);

% 可选：把被识别为“峰”的点画出来看看
plot(x(locs), A_final(locs), 'ko', 'MarkerSize',6, 'LineWidth',1.5);
