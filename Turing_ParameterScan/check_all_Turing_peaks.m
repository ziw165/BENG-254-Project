function check_all_Turing_peaks()
% 对所有 isTuring==true 的参数，一组一组跑 1D GM PDE，并统计峰数
% 需要：workspace 里已经有 params (N×4) 和 isTuring (N×1 logical)
% params 每行 = [aA, aI, rho, d]

    %----------------- 检查输入 -----------------%
    if ~evalin('base','exist(''params'',''var'')') || ...
       ~evalin('base','exist(''isTuring'',''var'')')
        error('请先在 base workspace 中生成 params 和 isTuring 再运行本函数。');
    end

    params   = evalin('base','params');
    isTuring = evalin('base','isTuring');

    idxT = find(isTuring);
    nT   = numel(idxT);

    fprintf('共有 %d 组参数满足图灵条件，现在逐一数峰...\n', nT);

    % 空间长度 L 要和你之前 PDE 模拟时保持一致
    L = 60;    % 你可以改成自己用的值

    peakCounts = zeros(nT,1);

    for ii = 1:nT
        p   = params(idxT(ii),:);
        aA  = p(1);
        aI  = p(2);
        rho = p(3);
        d   = p(4);

        fprintf('(%3d/%3d) aA=%.3g, aI=%.3g, rho=%.3g, d=%.3g ... ', ...
                ii, nT, aA, aI, rho, d);

        %----- 跑一次 1D PDE，取最终时刻的 A(x) -----%
        [x, Aend] = simulate_GM_1D(aA, aI, rho, d, L);

        %----- 用 findpeaks 数峰 -----%
        nPeaks = count_peaks(Aend, x);
        peakCounts(ii) = nPeaks;

        fprintf('peaks = %d\n', nPeaks);
    end

    % 把结果丢回 base workspace，方便你后续筛选
    assignin('base','idxTuring',idxT);
    assignin('base','peakCounts',peakCounts);

    % 例如：在 Command Window 里可以这样找 5 峰的参数：
    %   idx5 = idxTuring(peakCounts==5);
    %   params(idx5,:)

    fprintf('\n统计完成！结果已存入 base workspace: idxTuring, peakCounts。\n');
end


%================= 下面是局部函数 =================%

function [x, Aend] = simulate_GM_1D(aA, aI, rho, d, L)
% 用 pdepe 在 1D 上模拟无量纲 GM 模型
% dA/dt = aA*A^2/I - A + d2A/dx2
% dI/dt = aI*A^2   - rho*I + d*d2I/dx2

    Nx = 200;                   % 空间网格数
    x  = linspace(0, L, Nx);

    Tend  = 300;                % 终止时间
    Nt    = 100;                % 保存的时间点数
    tspan = linspace(0, Tend, Nt);

    % 均一稳态
    Astar = aA * rho / aI;
    Istar = aA^2 * rho / aI;

    noise_amp = 0.01;           % 初始噪声

    m = 0;  % 1D

    sol = pdepe(m, ...
        @(x,t,u,dudx) gm_pdefun(x,t,u,dudx,aA,aI,rho,d), ...
        @(x)          gm_icfun(x,Astar,Istar,noise_amp), ...
        @gm_bcfun, ...
        x, tspan);

    A = squeeze(sol(:,:,1));
    Aend = A(end,:);            % 最后一个时间点的 A(x)
end


function [c,f,s] = gm_pdefun(x,t,u,dudx,aA,aI,rho,d)
    A = u(1);  I = u(2);
    Ax = dudx(1); Ix = dudx(2);

    c = [1; 1];
    f = [Ax; d*Ix];             % D_A=1, D_I=d

    s1 = aA*A^2 / I - A;
    s2 = aI*A^2     - rho*I;
    s  = [s1; s2];
end


function u0 = gm_icfun(x,Astar,Istar,noise_amp)
    A0 = Astar * (1 + noise_amp * (2*rand-1));
    I0 = Istar * (1 + noise_amp * (2*rand-1));
    u0 = [A0; I0];
end


function [pl,ql,pr,qr] = gm_bcfun(xl,ul,xr,ur,t)
    % Neumann 边界：A_x = 0, I_x = 0
    pl = [0;0];
    ql = [1;1];
    pr = [0;0];
    qr = [1;1];
end


function nPeaks = count_peaks(Aend, x)
% 用 findpeaks 统计 A(x) 的峰数
    baseLevel = median(Aend);
    maxLevel  = max(Aend);
    minHeight = baseLevel + 0.2*(maxLevel - baseLevel);

    minDist = (x(end)-x(1))/20;  % 控制最小波长，可调整

    [pks, ~] = findpeaks(Aend, x, ...
                         'MinPeakHeight',   minHeight, ...
                         'MinPeakDistance', minDist);

    nPeaks = numel(pks);
end
