function nPeaks_pred = predict_peaks_GM(aA, ~, rho, d, L)
% 依据线性色散关系，预测在长度为 L 的 1D 区间上，模式的"峰数"
% 无量纲 GM:
% dA/dt = aA*A^2/I - A + d2A/dx2
% dI/dt = aI*A^2   - rho*I + d*d2I/dx2

    % 反应 Jacobian 元素（不含扩散）
    a  = 1;           % f_A
    b  = -1/aA;       % f_I
    c  = 2*aA*rho;    % g_A
    gv = -rho;        % g_I

    % 先确保本身就是 Turing 稳定的均一态（可选）
    trJ  = a + gv;
    detJ = a*gv - b*c;
    if ~(trJ < 0 && detJ > 0)
        nPeaks_pred = 0;
        return;
    end

    % 在一段连续的 k 上扫描最大实部特征值 lambda(k)
    kmax = 3.0;         % 最大波数（一般 2~3 足够）
    Nk   = 400;
    kvec = linspace(0, kmax, Nk);
    lam  = zeros(size(kvec));

    for i = 1:Nk
        k = kvec(i);
        a_k = a  - k^2;
        d_k = gv - d*k^2;
        tr_k  = a_k + d_k;
        det_k = a_k*d_k - b*c;
        disc  = tr_k^2 - 4*det_k;

        if disc >= 0
            lam(i) = 0.5*(tr_k + sqrt(disc));
        else
            lam(i) = real(0.5*(tr_k + 1i*sqrt(-disc)));
        end
    end

    [lam_max, idx] = max(lam);

    if lam_max <= 0
        % 没有 k>0 模式变不稳定，就不是 Turing
        nPeaks_pred = 0;
        return;
    end

    k_star = kvec(idx);         % 最不稳定波数

    % Neumann 边界下，约有 n ≈ k*L/(2π) 个峰
    nPeaks_pred = round( k_star * L / (2*pi) );
end
