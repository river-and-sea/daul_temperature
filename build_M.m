function M = build_M(kz_val, fr, param)
% build_M —— 8x10 边界条件矩阵（Mathematica 无量纲化简化版）
%   kz_val  轴向波数（标量）
%   fr      ki_omega 的输出（12 个频率相关量）
%   param   param_0 的静态量（b* 由 H_p/G_p 等现算；R、n 从 param 取）
%
% 列排布：[B_1, B_2, B_3, B_4, A_1, A_2, D_n, F_n, εni[1], εni[2]]
%   第1~8列 = 8 个力学未知量；第9、10列 = 井内入射项 εni[1]、εni[2] 的系数。
% 只输出 M（8x10）；需要方阵时调用方自己截 M(:,1:8)。

R = param.R;
n = param.n;

% ---- 从 fr / param 解包（变量名沿用原程序）----
ki_vec     = fr.ki_vec;
ks_val     = fr.ks_val;
kni_vec    = fr.kni_vec;
alpha_p_f  = fr.alpha1i;
alpha_p_ts = fr.alpha2i;
alpha_p_tf = fr.alpha3i;
alphas     = fr.alphas;
alpha_np   = fr.alphani;
kst        = fr.kst;
kft        = fr.kft;
knft       = fr.knft;

bH   = param.H_p/param.G_p;
bC   = param.C_p/param.G_p;
bPsi = param.Upsilon_p/param.G_p;
bO11 = param.Omega_11/param.G_p;
bO12 = param.Omega_12/param.G_p;
bO21 = param.Omega_21/param.G_p;
bO22 = param.Omega_22/param.G_p;
bMn  = param.Knf/param.G_p;
bOn  = param.Omega_2f_in/param.G_p;

hO11 = bO11*R^2; hO12 = bO12*R^2;hO21 = bO21*R^2;hO22 = bO22*R^2;hOn = bOn*R^2;

% ================================================================
% 第1步 & 第2步：计算径向波数与无量纲化 (对应 Mathematica 中的 khi, khni, khs 等)
% ================================================================
kri  = sqrt(kz_val^2 - ki_vec.^2);
knri = sqrt(kz_val^2 - kni_vec.^2);
krs  = sqrt(kz_val^2 - ks_val^2);

kappa_i  = kri  .* R;
kappa_ni = knri .* R;
kappa_s  = krs  * R;

khi  = ki_vec  .* R;  % 对应 Mathematica 中的 khi[j]
khni = kni_vec .* R;  % 对应 Mathematica 中的 khni[i]
khs  = ks_val  * R;   % 对应 Mathematica 中的 khs
khz  = kz_val  * R;   % 对应 Mathematica 中的 khz

% ================================================================
% 第3步：计算 Bessel 函数比值（此部分处理不同分支，严格保留）
%   n_zhi / n_zhi_s：=1 时在 K_n 比值里加 -1i*pi*I_n 修正项（处理复数宗量跨割线）。
%   ★★ 取值已实测定案（新全波_check\nz_test.m、nz_probe_run.m，2026-10-10）：
%     本算例 n=1/pc=10MPa 全频带内 Re κ 从不为负（S 支 min 0.0187/168100 次、
%     井外 4 支 min 0.0093/672400 次，Re<0 计数皆 0）⇒ 不跨割线 ⇒ 修正项不该加；
%     实测改 n_zhi_s=1 会把 max|p| 压到 8.8e-4、与刘程序 corr 掉到 0.0288；故必须为 0（与 build_M_jvre.m 一致），旧注"必须为 1"已被推翻。
% ================================================================
n_zhi    = [0, 0, 0, 0];   % 井外 4 支 P 波：0 = 纯 K_n（原程序同）
n_zhi_s  = 0;              % S 波支：0 = 纯 K_n（实测 Re κ 不跨割线；改成 1 波形崩）
n_zhi_ni = [0, 0];         % 井内 2 支热波（第 9、10 列的 K_n 侧）：0 = 纯 K_n

CalK_i = kappa_i(:) .* ...
    (besselk(n+1, kappa_i(:), 1) - 1i*pi*n_zhi(:).*besseli(n+1, kappa_i(:), 1)) ./ ...
    (besselk(n,   kappa_i(:), 1) - 1i*pi*n_zhi(:).*besseli(n,   kappa_i(:), 1));

CalK_s = kappa_s * ...
    (besselk(n+1, kappa_s, 1) - 1i*pi*n_zhi_s*besseli(n+1, kappa_s, 1)) ./ ...
    (besselk(n,   kappa_s, 1) - 1i*pi*n_zhi_s*besseli(n,   kappa_s, 1));

CalI_ni = kappa_ni(:) .* besseli(n+1, kappa_ni(:), 1) ./ besseli(n, kappa_ni(:), 1);

% ================================================================
% 第4步：定义 Mathematica 替换规则对应的中间变量
% ================================================================
KK  = n - CalK_i;       % 对应 Mathematica 中的 KK[j]
II  = n + CalI_ni;      % 对应 Mathematica 中的 II[i]
KKs = n - CalK_s;       % 对应 Mathematica 中的 KKs

% ================================================================
% 第5步：直接按 Mathematica 矩阵结构向量化赋值 (整体已除以 R^2)
% 列排布：[B_1, B_2, B_3, B_4,  A_1, A_2,  D_n, F_n,  εni[1], εni[2]]
% ================================================================
M = zeros(8, 10);

% 行1：法向应力连续
M(1, 1:4) = 2*n^2 - 2*KK - khi.^2 .* (-2 + bH + bC .* alpha_p_f) ...
            - hO11 .* alpha_p_ts - hO12 .* alpha_p_tf + 2 * kappa_i.^2;
M(1, 5:6) = bMn * khni.^2 + hOn .* alpha_np;
M(1, 7)   = 2 * (-1 + KKs) * n;
M(1, 8)   = 2 * (-KKs + n^2 + kappa_s^2);

% 行2：孔隙流体压力连续
M(2, 1:4) = -khi.^2 .* (bC + bPsi .* alpha_p_f) ...
            - hO21 .* alpha_p_ts - hO22 .* alpha_p_tf;
M(2, 5:6) = bMn * khni.^2 + hOn .* alpha_np;
% M(2, 7:8) 已经是 0

% 行3：切向应力连续 (τ_{rθ})
M(3, 1:4) = n - n .* KK;
% M(3, 5:6) 已经是 0
M(3, 7)   = KKs - n^2 - kappa_s^2 / 2;
M(3, 8)   = n - KKs * n;

% 行4：切向应力连续 (τ_{rz})
M(4, 1:4) = KK;
% M(4, 5:6) 已经是 0
M(4, 7)   = n / 2;
M(4, 8)   = KKs * (1- (khs/khz)^2/2);

% 行5：径向位移连续 (u_r，此处 tiaojian_i=1 已代入)
M(5, 1:4) = -KK .* (1 + alpha_p_f);
M(5, 5:6) = II;
M(5, 7)   = -n * (1 + alphas);
M(5, 8)   = -KKs * (1 + alphas);

% 行6：流体温度连续 (T_f)
M(6, 1:4) = -alpha_p_tf;
M(6, 5:6) = alpha_np;

% 行7：固体温度连续 (T_s)
M(7, 1:4) = -alpha_p_ts;
M(7, 5:6) = alpha_np;

% 行8：热流连续 (q_r)
M(8, 1:4) = KK .* (kst .* alpha_p_ts + kft .* alpha_p_tf);
M(8, 5:6) = -knft .* II .* alpha_np;
% M(8, 7:8) 已经是 0

% ================================================================
% 第6步：第9、10列 —— 井内入射项 εni[1]、εni[2] 的系数
% ================================================================
CalK_ni = kappa_ni(:) .* ...
    (besselk(n+1, kappa_ni(:), 1) - 1i*pi*n_zhi_ni(:).*besseli(n+1, kappa_ni(:), 1)) ./ ...
    (besselk(n,   kappa_ni(:), 1) - 1i*pi*n_zhi_ni(:).*besseli(n,   kappa_ni(:), 1));
KKn = n - CalK_ni;                          % = n - κn K_{n+1}(κn)/K_n(κn)，与 II = n + În[i] 对偶

M(1, 9:10) = bMn*khni.^2 + hOn.*alpha_np;   % 行1：与第5、6列同式
M(2, 9:10) = bMn*khni.^2 + hOn.*alpha_np;   % 行2：与第5、6列同式
M(5, 9:10) = KKn;                           % 行5：u_r 的 K_n 侧组合
M(8, 9:10) = -knft .* KKn .* alpha_np;      % 行8：热流的 K_n 侧组合
% 行3、4、6、7 恒为 0：第 9、10 列只出现在第 1、2、5、8 行。
% （激发留数那一步会用到这两列，det 只用前 8 列。）

% ================================================================
% 弹性退化
% ================================================================
% if n == 0
%     M1 = M([1,3:8], [1,3:8]);
% else
%     M1 = M;
% end
% 
% M_r = real(M1);
% M_i = imag(M1);5
% M_x = abs(M1);

end