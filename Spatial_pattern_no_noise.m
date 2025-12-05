%% Part C: Turing-like digit pattern in 1D using activator-inhibitor RD

clear; clc; close all;

%% Spatial domain
N = 600;          % number of spatial points
L = 2.3;          % limb length
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
A0 = 0.1*(1 + 0.01*ones(N,1));
I0 = 0.1* (1 + 0.01*ones(N,1));

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
