function out = axis_fullwave(n, pc_MPa, n_rec, theta0)
% axis_fullwave —— 井轴上的"广义声压"全波列（张碧星等 1995，式 (9)）
%
%   只做两件事：逐 (ω, kz) 求出式 (9) 第二项的被积函数 p_r，然后做两次 FFT
%   得到时域总场 P(z,t) = 直达项 + 反射项。
%
%   参数/网格都是本文件内的常量（原 fw_*），介质与几何来自 param_0。
%   井内两支：M8_solve 给出 An（解级、按 I_n 归一），只有第 1 支（声支）参与
%   "广义声压"；热支的贡献在温度场里单独算。
%
%   ── 附注（与旧 zbx_axis_pressure 的差别）──
%   旧程序的 src 与 An 都带 exp(-|Re κ|)；现在 M8_solve 里都不带，
%   故 p_r 与该旧程序的数值会有整体口径差，但 (ω,k) 上的相对形状一致。

if nargin < 2, pc_MPa = 0; end
if nargin < 3, n_rec  = 3;  end
if nargin < 4, theta0 = 0;  end

% 近奇异点警告静默
w1 = warning('off','MATLAB:singularMatrix');
w2 = warning('off','MATLAB:illConditionedMatrix');
w3 = warning('off','MATLAB:nearlySingularMatrix');
restoreWarn = onCleanup(@() warning([w1 w2 w3]));   %#ok<NASGU>

%% ============ 0. 参数 ============
pc = pc_MPa*1e6;
pf = 0;
R  = 0.1;
param = param_0(pc, pf, theta0, R, n);

G = grid_setup(n, n_rec);
r0    = param.a0;
s0    = 1e-6;
eps_n = G.eps_n;
angf  = 1;                                   % cos(n(θ-θ0))，θ=θ0 时 = 1

%% ============ 1. 逐 (ω, kz) 求式 (9) 的被积函数 + 直达项 ============
p_r  = zeros(G.n_w, G.n_k);
A_ax = zeros(G.n_w, G.n_k);
P_d  = zeros(G.n_w, numel(G.z_rec));
nbad = 0;  kmax = 0;
chengzi = angf * eps_n * (r0/4)^n / (4*pi*factorial(n));

for iw = 1:G.n_w
    fr   = ki_omega(G.om_c(iw), param, []);  % 复频下的解
    k_in = fr.kni_vec(:);                    % 2x1 井内两支轴向波数（列向量，直达项要做外层积）

    % ---- 直达项（式 (9) 第一项，两支求和）----
    Dn = zeros(1, numel(G.z_rec));
    for m = 0:n
        Dn = Dn + nchoosek(n,m)*factorial(2*n-m) ...
            .* sum(exp(1i*(k_in*G.z_rec)) .* ((-2i*k_in).^m) ./ (G.z_rec.^(2*n-m+1)), 1);
    end
    P_d(iw,:) = chengzi .* Dn;

    % ---- 反射项的被积函数 ----
    for ik = 1:G.n_k
        kz  = G.k_z(ik);
        knr = sqrt(kz.^2 - k_in.^2);         % 2x1 井内径向波数
        kap = knr*R;

        [~, ~, An_o] = M8_solve(fr, param, kz);      % 只取解级系数（不算 Dp）
        An = An_o.' ./ besseli(n, kap, 1);           % 2x1，按 I_n 归一

        p_r(iw,ik)  = chengzi * sum((knr.^(2*n)) .* An);
        A_ax(iw,ik) = An(1);
        if ~isfinite(p_r(iw,ik)), p_r(iw,ik) = 0; nbad = nbad + 1; end
        kmax = max(kmax, max(abs(kap)));
    end
end

%% ============ 3. 双重 FFT -> 时域 ============
P_r_w = integrate_kz(p_r, G) .* G.F(:);       % 反射项：波数变换后乘源谱
P_d_w = P_d .* G.F(:) / G.dt;                 % 直达项：乘源谱及采样时间间隔

T_r = integrate_omega(P_r_w, G);
T_d = integrate_omega(P_d_w, G);
decay = exp(-G.w0*G.i_t*G.dt).';
P_r_t = real(T_r(G.i_t, :) .* decay);
P_d_t = real(T_d(G.i_t, :) .* decay);
P = (P_r_t + P_d_t).';

out = struct('t', G.t_ms, 'P', P, 'z', G.z_rec, 'p_r', p_r, 'A_ax', A_ax, ...
             'Pr_t', P_r_t.', 'Pd_t', P_d_t.', ...
             'n', n, 'pc', pc, 'fc', G.f_src, ...
             'kz', G.k_z, 'freq', G.freq, 'nbad', nbad, 'kapmax', kmax);

fprintf('n=%d pc=%g MPa  z=%.4f..%.4f m  max|P_r|=%.4g  max|P_d|=%.4g  max|P|=%.4g  非有限点=%d\n', ...
    n, pc_MPa, G.z_rec(1), G.z_rec(end), max(abs(out.Pr_t(:))), max(abs(out.Pd_t(:))), ...
    max(abs(P(:))), nbad);
end


% ======================================================================
function G = grid_setup(n_src, n_rec)
% 频率/波数/时间三套栅格 + 源谱（常量沿用原 fw_*）
fw_tc      = 0.5e-3;
fw_fc_list = [6000 5000 6000];
fw_n_list  = [0 1 2];
fw_dt      = 10e-6;
fw_n2      = 2049;
fw_n1      = 4097;
fw_dz      = 0.5;
fw_ddz0    = 0.0208;
fw_zp      = 2;
fw_f0i_ratio = 1/100;

G.f_src   = fw_fc_list(fw_n_list == n_src);
G.eps_n   = 1 + (n_src > 0);
tc        = fw_tc;
G.f0      = -G.f_src*fw_f0i_ratio;
G.w0      = 2*pi*G.f0;

G.dt      = fw_dt;   n1 = fw_n1;   n2 = fw_n2;
n_step    = fix(fw_dz/fw_ddz0);
G.d_z     = fw_dz/n_step;
G.n_step  = n_step;
G.dk      = 2*pi/(G.d_z*(n1-1));
m0        = round(fw_zp/G.d_z);
G.m0      = m0;
G.m_rec   = 1 + m0 + (0:(n_rec-1))*n_step;
G.z_rec   = G.m_rec*G.d_z;
G.n_k     = fix(n1/2) + 1;
G.k_z     = ((1:G.n_k)-1)*G.dk + 1e-36;
G.n_t     = fix(0.5*(fix(n2/2)+1));
G.i_t     = 1:G.n_t;
G.t_ms    = (G.i_t)*G.dt*1000 - G.dt*1000;

G.dfp     = 1/(G.dt*(n2-1));
G.iw_f    = (round((G.f_src-2/tc)/G.dfp)+1) : round((G.f_src+2/tc)/G.dfp);
G.freq    = G.iw_f*G.dfp;
G.freq_c  = G.freq - 1i*G.f0;
G.om_c    = 2*pi*G.freq_c;
G.F       = source_spec(G.freq_c, G.f_src, tc);
G.n_w     = numel(G.iw_f);
G.n1 = n1;  G.n2 = n2;
end


% ======================================================================
function Zz = integrate_kz(Z, G)
% 对矩阵列方向做波数域变换并选取接收道
Zk = zeros(size(Z,1), G.n1, 'like', Z);
Zk(:, 1:G.n_k)          = Z;
Zk(:, G.n1+1-(2:G.n_k)) = Z(:, 2:G.n_k);      % 负半轴：复制式镜像
rec_z = fft(Zk(:, 1:G.n1-1).').';             % kz -> z
Zz = rec_z(:, G.m_rec);
end


% ======================================================================
function Zt = integrate_omega(Z, G)
% 补负频共轭并执行频域到时间域的 FFT
Xf = zeros(G.n2-1, size(Z,2), 'like', Z);
Xf(G.iw_f, :)        = Z;                     % 正频
Xf(G.n2-G.iw_f+1, :) = conj(Z);               % 负频共轭镜像
Zt = fft(Xf, [], 1);                          % omega -> t
end
