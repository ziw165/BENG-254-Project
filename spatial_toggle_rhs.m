function dydt = spatial_toggle_rhs(t, y, N, dx, DA, DI, alphaA, alphaI, m, n)
% Reaction–diffusion Sox9–Wnt toggle model on a 1D grid
% y = [A1..AN, I1..IN]

A = y(1:N);
I = y(N+1:end);

dA = zeros(N,1);
dI = zeros(N,1);

% 内点：i = 2..N-1，中心差分
for i = 2:N-1
    lapA = (A(i+1) - 2*A(i) + A(i-1))/dx^2;
    lapI = (I(i+1) - 2*I(i) + I(i-1))/dx^2;
    
    dA(i) = alphaA/(1 + I(i)^m) - A(i) + DA*lapA;
    dI(i) = alphaI/(1 + A(i)^n) - I(i) + DI*lapI;
end

% 左边界：Neumann（zero-flux）近似
lapA1 = (A(2) - A(1))/dx^2;
lapI1 = (I(2) - I(1))/dx^2;
dA(1) = alphaA/(1 + I(1)^m) - A(1) + DA*lapA1;
dI(1) = alphaI/(1 + A(1)^n) - I(1) + DI*lapI1;

% 右边界：Neumann
lapAN = (A(N-1) - A(N))/dx^2;
lapIN = (I(N-1) - I(N))/dx^2;
dA(N) = alphaA/(1 + I(N)^m) - A(N) + DA*lapAN;
dI(N) = alphaI/(1 + A(N)^n) - I(N) + DI*lapIN;

dydt = [dA; dI];
end
