function flag = isTuringGM(aA, aI, rho, d, L, nMaxK)
% flag = true 说明该参数组存在 k>0 的空间模态不稳定（图灵斑纹）
% 无量纲 GM 模型：
%   dA/dt = aA*A^2/I - A + d2A/dx2
%   dI/dt = aI*A^2   - rho*I + d*d2I/dx2

%---------- 1. 均一稳态 ----------
Astar = aA * rho / aI;
Istar = aA^2 * rho / aI;
% 这里没有显式用到，但如果以后要做非线性模拟，可以直接复用

%---------- 2. 反应 Jacobian (不含扩散) ----------
% f(A,I) = aA*A^2/I - A
% g(A,I) = aI*A^2   - rho*I
% 在稳态处的偏导，利用稳态关系可大幅简化：
a  = 1;           % f_A
b  = -1/aA;       % f_I
c  = 2*aA*rho;    % g_A
dd = -rho;        % g_I   (避免与扩散比 d 重名)

%---------- 3. k = 0 模式必须稳定（无扩散稳定） ----------
trJ  = a + dd;
detJ = a*dd - b*c;     % det(J)

if trJ >= 0 || detJ <= 0
    flag = false;      % 无扩散都不稳定，就不是"图灵斑纹"
    return;
end

%---------- 4. 扫 k>0 的空间模态，寻找扩散诱导的不稳定 ----------
flag = false;
for n = 1:nMaxK
    k = n*pi/L;        % Neumann 边界的波数

    a_k = a  - (k^2);      % 因为 D_A = 1
    d_k = dd - d*(k^2);    % 抑制剂扩散系数 = d

    tr_k  = a_k + d_k;
    det_k = a_k*d_k - b*c;

    % 最大特征值
    disc = tr_k^2 - 4*det_k;
    if disc >= 0
        lambda_max = 0.5*(tr_k + sqrt(disc));
    else
        lambda_max = 0.5*(tr_k + 1i*sqrt(-disc));
    end

    if real(lambda_max) > 0
        flag = true;       % k>0 模式不稳定 => 图灵不稳定
        return;
    end
end
end
