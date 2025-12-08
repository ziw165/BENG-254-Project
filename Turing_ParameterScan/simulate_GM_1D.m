function [x, Aend] = simulate_GM_1D(aA, aI, rho, d, L)
% 用 pdepe 在 1D 上模拟无量纲 GM 模型

Nx = 200;
x = linspace(0,L,Nx);

Tend = 300;
Nt = 100;
tspan = linspace(0,Tend,Nt);

Astar = aA * rho / aI;
Istar = aA^2 * rho / aI;

noise = 0.01;

m = 0;
sol = pdepe(m, ...
    @(x,t,u,dudx) gm_pde(x,t,u,dudx,aA,aI,rho,d), ...
    @(x) [Astar*(1+noise*(2*rand-1));
          Istar*(1+noise*(2*rand-1))], ...
    @gm_bc, ...
    x, tspan);

A = squeeze(sol(:,:,1));
Aend = A(end,:);
end


function [c,f,s] = gm_pde(x,t,u,dudx,aA,aI,rho,d)
A = u(1); I = u(2);
Ax = dudx(1); Ix = dudx(2);

c = [1;1];
f = [Ax; d*Ix];
s = [aA*A^2/I - A;
     aI*A^2   - rho*I];
end


function [pl,ql,pr,qr] = gm_bc(xl,ul,xr,ur,t)
pl = [0;0]; ql = [1;1];
pr = [0;0]; qr = [1;1];
end
