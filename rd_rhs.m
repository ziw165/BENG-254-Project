function dYdt = rd_rhs(~,Y,N,dx,DA,DI,alphaA,alphaI,muA,muI)

    % 拆成 A, I
    A = Y(1:N);
    I = Y(N+1:end);

    % ---- 1. 防止出现负浓度 ----
    A(A < 0) = 0;
    I(I < 0) = 0;

    % ---- 2. 给 I 一个下限，避免 A.^2 ./ I 爆炸 ----
    I_floor = 1e-3;              % 可以试 1e-4 ~ 1e-3
    I_safe  = I;
    I_safe(I_safe < I_floor) = I_floor;

    % ---- 3. 计算 Laplacian，带 Neumann BC（向量化写法）----
    LapA = zeros(N,1);
    LapI = zeros(N,1);

    % 中间点
    LapA(2:N-1) = (A(3:N)   - 2*A(2:N-1)   + A(1:N-2))   / dx^2;
    LapI(2:N-1) = (I(3:N)   - 2*I(2:N-1)   + I(1:N-2))   / dx^2;

    % 边界：零通量（复制邻点）
    LapA(1) = LapA(2);
    LapA(N) = LapA(N-1);
    LapI(1) = LapI(2);
    LapI(N) = LapI(N-1);

    % ---- 4. 反应项 + 扩散 ----
    dA = alphaA .* (A.^2 ./ I_safe) - muA .* A + DA .* LapA;
    dI = alphaI .* (A.^2)           - muI .* I + DI .* LapI;

    dYdt = [dA; dI];
end
