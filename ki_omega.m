function fr = ki_omega(omega, param, prev_lambda)
% ki_omega —— 单个频点：求根 → 沿支追踪 → 全部随 ω 变化的量
%
%   omega        角频率 [rad/s]
%   param        param_0(pc, pf, theta0) 的结果
%   prev_lambda  上一频点的 4 个根（首点传 []，此时按实部排序）
%
% 输出：
%   12 个量，全部直接进 build_M_wulianggang / build_M_jvre / compute_det
%     alpha_s_out  ki_vec  ks_val  kni_vec  alpha1i  alpha2i  alpha3i  alphas  alphani
%     kst  kft  knft
%   lambda_current  本频点的 4 个根（配好支序），下一频点作为 prev_lambda
%
% 符号句柄：首次调用时由 param_sym 建好（persistent），之后全是数值求值。

persistent f_poly f_jie f_det
if isempty(f_poly)
    [C_sym, jie_sym, det_sym, vars] = param_sym();
    f_poly = matlabFunction(C_sym,  'Vars', vars(2:end));
    f_jie  = matlabFunction(jie_sym, 'Vars', vars);
    f_det  = matlabFunction(det_sym, 'Vars', vars);
end

den = param.den;  den_f = param.den_f;  T0 = param.T0;

den_tilde = 1i*param.visc/omega/param.k0 + den_f*param.morsE/param.poros;
F_xz      = param.K/T0/(-1i*omega);

% ---- 16 个本构量（各自具名）----
H_p       = param.H_p;
C_p       = param.C_p;
Upsilon_p = param.Upsilon_p;
K_s_tilde = param.k_s/T0/(1i/omega + param.time0_s);
K_f_tilde = param.k_f/T0/(1i/omega + param.time0_f);
Omega_11  = param.Omega_11;
Omega_12  = param.Omega_12;
Omega_21  = param.Omega_21;
Omega_22  = param.Omega_22;
F_d11     = param.F_p11 + F_xz;
F_d12     = param.F_p12 - F_xz;
F_d21     = param.F_p21 - F_xz;
F_d22     = param.F_p22 + F_xz;

cons = {H_p, C_p, Upsilon_p, K_s_tilde, K_f_tilde, ...
        Omega_11, Omega_12, Omega_21, Omega_22, ...
        den, den_f, den_tilde, F_d11, F_d12, F_d21, F_d22};

% ---- 求根 + 沿支追踪 ----
S_squared = track_lambda(roots(f_poly(cons{:})), prev_lambda);

% ---- 井外 4 支 ----
ki_vec = omega*sqrt(S_squared);
ks_val = omega*sqrt((den - den_f^2/den_tilde)/param.G_p);

% ---- 井内 2 支：二次方程解析根 ----
K_f_tilde_in = param.k_f_in/T0/(1i/omega + param.time0_f_in);
b_2jie = param.den_f_in*K_f_tilde_in + param.Knf*param.F_p22_in + param.Omega_2f_in^2;
disc   = sqrt(b_2jie^2 - 4*param.Knf*K_f_tilde_in*param.den_f_in*param.F_p22_in);
S_squared_in = [(b_2jie - disc)/(2*param.Knf*K_f_tilde_in); ...
                (b_2jie + disc)/(2*param.Knf*K_f_tilde_in)];
kni_vec = omega*sqrt(S_squared_in);

% ---- 井外两相系数 ----
alpha_blk = f_jie(S_squared(:).', cons{:}) ./ ...
            f_det(S_squared(:).', cons{:});

alpha1i    = alpha_blk(1,:).';                  % 不乘 ω²
alpha2i    = omega^2*alpha_blk(2,:).';
alpha3i    = omega^2*alpha_blk(3,:).';
alpha_s_out = -den_f/den_tilde;
alphas     = alpha_s_out;
alphani    = (den_f - param.Knf*S_squared_in)*omega^2/param.Omega_2f_in;

% ---- 热导率比 ----
kst  = param.k_s   /(1 - 1i*param.time0_s   *omega);
kft  = param.k_f   /(1 - 1i*param.time0_f   *omega);
knft = param.k_f_in/(1 - 1i*param.time0_f_in*omega);

fr = struct('alpha_s_out', alpha_s_out, 'ki_vec', ki_vec, 'ks_val', ks_val, ...
            'kni_vec', kni_vec, 'alpha1i', alpha1i, 'alpha2i', alpha2i, ...
            'alpha3i', alpha3i, 'alphas', alphas, 'alphani', alphani, ...
            'kst', kst, 'kft', kft, 'knft', knft, ...
            'lambda_current', S_squared);
end
