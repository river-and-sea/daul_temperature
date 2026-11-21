function param = param_0(pc, pf, dT, n, r, material)
% param_0 —— 不随频率变化的全部参数（算例量 + 材质 + 派生量）
%
%   param = param_0()                          % 全默认（pc=65 MPa, pf=5 MPa, dT=0, n=0, r=R，赵组材质）
%   param = param_0(pc, pf, dT, n, r)          % 位置参数（可空取默认）
%   param = param_0(pc, pf, dT, n, r, 'chen')  % 第 6 个可选：材质预设名
%
%   输入
%     pc  围压 [Pa]             默认 65e6
%     pf  孔压 [Pa]             默认 5e6
%     dT  温度变化 [K]          默认 0
%     n   方位序数（Bessel 阶） 默认 0
%     r   接收器径向位置 [m]    默认 R（井壁）
%     material  材质预设名      默认 'zhao'（可选 'chen' / 'hou'）
%
%   材质结构体（mat_preset 的返回）
%     mat.elastic  孔隙弹性 + VTI 各向异性
%     mat.thermal  热学与耗散
%
%   输出
%     param.case   算例量 struct('pc','pf','dT','n','r')
%     其余字段     材质派生量（A3333 / G_p / ... / R / r0 / s0 / th_rec / th_src）

% ---------------- 算例量 ----------------
R  = 0.1;                                         % 井半径 [m]（固定，不作输入）
if nargin < 1 || isempty(pc), pc = 0e6; end
if nargin < 2 || isempty(pf), pf = 0e6;  end
if nargin < 3 || isempty(dT), dT = 10;    end
if nargin < 4 || isempty(n),  n  = 0;    end
if nargin < 5 || isempty(r),  r  = R;    end
if nargin < 6 || isempty(material), material = 'zhao'; end

% ---------------- 材质（弹性 + 热学）----------------
mat = mat_preset(material);
e   = mat.elastic;
th  = mat.thermal;

% ---------------- 流体黏度 ----------------
% ★ 2026-11-22 用户要求：黏度**固定**，不随温度变化（visc = visc0）
%   原先的 Vogel 型温度修正 μ(T)=μ0·10^[247.8/(T-140) - 247.8/(T0-140)] 已停用，
%   需要恢复时把下面两行注释打开、并删掉 visc = e.visc0 那一行即可。
Tp = th.T0 + dT;
% if Tp <= 160
%     error('param_0:dT', '温度 %.1f K 过低（Vogel 公式要求 T > 140 K）', Tp);
% end
% visc = e.visc0 * 10^(247.8/(Tp - 140) - 247.8/(th.T0 - 140));
visc = e.visc0;

% ---------------- 源径与方位角（几何）----------------
r0 = 0.01;                    % 源分布半径 [m]
s0 = 1e-9;                    % 源强标定（原 fw_src_amp）
th_rec = 0;                   % 接收方位角 [rad]
th_src = 0;                   % 源方位角 [rad]

%% §1 井外介质派生量
% —— Biot 孔弹性模量 ——
Alpha    = 1 - e.Kb/e.Ks;                          % Biot 有效应力系数 α（0 < α < 1）
M        = 1/((Alpha - e.poros)/e.Ks + e.poros/e.Kf);  % Biot 模量 M [Pa]
lambda_b = e.Kb - 2/3*e.mu_b;                      % 干岩石 Lamé 第一参数 λ_b [Pa]
lambda_c = lambda_b + M*Alpha^2;                   % 不排水/封闭 λ_c [Pa]
Kc       = e.Kb + M*Alpha^2;                       % 封闭体积模量 Kc [Pa]

% —— 相密度 ——
den_s_prime = e.den_s*(1 - e.poros);               % (1-φ)ρs
den_f_prime = e.den_f*e.poros;                     % φρf
den = 2270;      % 总密度：原程序口径固定 2270

% —— 热惯性系数（= 等效密度 × 比热容）——
F11 = den_s_prime*th.C_E_s;
F22 = den_f_prime*th.C_E_f;
F12 = 0;  F21 = 0;                                 % 两相热惯性不交叉

% —— 热-孔弹性耦合系数（温度变化 → 应力 / 孔压）——
Omega_1s = e.poros*M*(1-e.poros)/(1-Alpha)*(Kc/e.poros/M-Alpha)*th.alpha_s;
Omega_1f = e.poros*M*th.alpha_f*Alpha;
Omega_2s = e.poros*M*(1-e.poros)/(1-Alpha)*(Alpha/e.poros-1)*th.alpha_s;
Omega_2f = e.poros*M*th.alpha_f;

% —— 温度变化与热惯性（除以当前温度 Tp=T0+dT，直接进 8×8 矩阵）——
theta_s_s = dT;  theta_f_s = dT;                   % 固/流两相取同一温度变化 ΔT
F_p11 = F11/Tp;  F_p12 = F12/Tp;
F_p21 = F21/Tp;  F_p22 = F22/Tp;

% —— 相热导率 ——
k_f = e.poros*th.k0_f;
k_s = (1 - e.poros)*th.k0_s;

% —— 井内流体（井内 2 支模式用；口径：Knf = Kf，不引入 Q 衰减）——
Knf = e.Kf;  k_f_in = th.k0_f;
% 井内静态压缩 w'_内（(70) 式，平衡井 P_内 = pf），只用于迹因子修正
w_in = -(pf - th.alpha_f*Knf*dT)/(3*Knf);
% Ω_内2 的迹因子修正：Ω'_内2 = Ω_内2 (1 + 3 w'_内)（对应表末行）
Omega_2f_in = th.alpha_f*Knf*(1 + 3*w_in);
den_f_in = e.den_f;
F_p22_in = den_f_in*th.C_E_f/Tp;
time0_f_in = th.time0_f;

%% §2 热致有效压力（静态温升经 Omega 折算为等效压力，与围压/孔压线性叠加）
%  pc_t —— 热致有效围压：外压 - 温升项（固相侧 Omega_1s、流体侧 Omega_1f）
%  pf_t —— 热致有效孔压：孔压 - 温升项（固相侧 Omega_2s、流体侧 Omega_2f）
pc_t =  pc - Omega_1s*theta_s_s - Omega_1f*theta_f_s;
pf_t =  pf - Omega_2s*theta_s_s - Omega_2f*theta_f_s;

%% §3 静态形变
e_s    = (-pc_t + Alpha*pf_t)/e.Kb;               % 骨架体积应变
zeta_s = (-Alpha*pc_t + (lambda_c/M)*pf_t)/e.Kb;  % 孔隙体积应变增量
e33    = e_s/3;
w33    = -zeta_s/3;
AlphaM = Alpha*M;

%% §4 四个模量（干岩石模量 + 三阶常数修正）
A3333 = lambda_c + 2*e.mu_b ...
        + (lambda_c + e.nu1 + 2*e.nu2)*e_s ...
        + 2*(lambda_c + 3*e.mu_b + 2*e.nu2 + 4*e.nu3)*e33 ...
        - zeta_s*(AlphaM - 2*e.gama2);

G_p = e.mu_b + (lambda_c + e.nu2 + (4/3)*(e.mu_b + e.nu3))*e_s ...
      - (AlphaM + e.gama/2)*zeta_s;

C_p = AlphaM + (AlphaM + e.gama)*e33 ...
      - (2*e.gama2 + e.gama)*e_s ...
      - zeta_s*(AlphaM + 2*e.gama3);

Upsilon = M + e_s*(2*e.gama3 + AlphaM) ...
          - zeta_s*3*(M - 2*e.gama1);

% —— 再扣掉热应力项 ——
H_p       = A3333   - (Omega_1f*theta_f_s + Omega_1s*theta_s_s);
G_p       = G_p     - (Omega_1f*theta_f_s + Omega_1s*theta_s_s);
Upsilon_p = Upsilon - (Omega_2f*theta_f_s + Omega_2s*theta_s_s);

% —— 动态耦合系数（带 (1+形变) 因子）——
Omega_11 = Omega_1s*(1 + e33);  Omega_12 = Omega_1f*(1 + e33);
Omega_21 = Omega_2s*(1 + 3*w33);  Omega_22 = Omega_2f*(1 + 3*w33);

%% 打包（算例量单独进 param.case，其余字段名不变）
param = struct( ...
    'case',      struct('pc',pc, 'pf',pf, 'dT',dT, 'n',n, 'r',r), ...
    'A3333',     A3333,      'G_p',       G_p,       ...
    'C_p',       C_p,        'Upsilon',   Upsilon,   ...
    'H_p',       H_p,        'Upsilon_p', Upsilon_p, ...
    'Omega_11',  Omega_11,   'Omega_12',  Omega_12,  ...
    'Omega_21',  Omega_21,   'Omega_22',  Omega_22,  ...
    'F_p11',     F_p11,      'F_p12',     F_p12,     ...
    'F_p21',     F_p21,      'F_p22',     F_p22,     ...
    'k_f',       k_f,        'k_s',       k_s,       ...
    'time0_s',   th.time0_s, 'time0_f',   th.time0_f, ...
    'T0',        th.T0,      'Tp',        Tp,        'K',         th.K,      ...
    'den',       den,        'den_f',     e.den_f,   ...
    'poros',     e.poros,    'k0',        e.k0,      ...
    'morsE',     e.morsE,    'visc',      visc,      ...
    'k_f_in',    k_f_in,     'Knf',       Knf,       ...
    'Omega_2f_in', Omega_2f_in, 'den_f_in', den_f_in, ...
    'F_p22_in',  F_p22_in,   'time0_f_in', time0_f_in, ...
    'R',         R,          'r0',        r0,        ...
    's0',        s0,         'th_rec',    th_rec,    ...
    'th_src',    th_src);
end

% ======================================================================
function mat = mat_preset(name)
% mat_preset —— 材质预设（每组自带「弹性 mat.elastic + 热学 mat.thermal」）
switch lower(name)
    case 'zhao'      % 赵志强（原默认）
        mat.elastic = struct( ...
            'k0',    43.5e-15, 'morsE', 2.8481, 'visc0', 8.5e-4, ...
            'den_s', 2581.57,  'poros', 0.197,  'den_f', 1000, ...
            'Kf',    2.25e9,   'Ks',    37e9,   ...
            'mu_b',  5.4e9,    'Kb',    13.2e9, ...
            'nu1',   3578e9,   'nu2',  -1310e9, 'nu3',  -345e9, ...
            'gama', -1150e9,   'gama1',  1.93e9, ...
            'gama2',  220e9,   'gama3', -110e9);
        mat.thermal = struct( ...
            'T0',      300,    'k0_s',   7.69,   'k0_f',   0.598, ...
            'time0_s', 1e-11,  'time0_f', 1e-9, ...
            'C_E_s',   740,    'C_E_f',  4157, ...
            'K',       1e-7,   'alpha_s', 3.45e-5, 'alpha_f', 2.07e-4);
        % ==== 全部为实测/原表值（20 degC 基准），出处见 06b 附录 A.5.5 ====
        %   k0_s  = 7.69     Horai (1971) JGR 76(5) 1278-1308：多晶石英聚集体实测
        %   k0_f  = 0.598    NIST SRD 69（IAPWS R15-11），纯水 20 degC 饱和线
        %   C_E_s = 740      Waples & Waples (2004) Nat.Resour.Res. 13(2)：石英量热值 20 degC
        %   C_E_f = 4157     NIST SRD 69（IAPWS-95），纯水 c_v 20 degC 饱和线
        %   alpha_s = 3.45e-5 Kosinski et al. (1992) EFTF Table 1（HTMIAC 实测推荐集），体膨胀
        %   alpha_f = 2.07e-4 Engineering ToolBox 体膨胀原表值，水 20 degC
        % ↑ 20 degC 基准的文献值（出处见 06b 附录 A.5.4）
        %   k0_s=7.71  石英单晶体平均（Voigt），与多晶实测 7.69 差 0.3%
        %   k0_f=0.598 纯水 20 degC（NIST）
        %   C_E_s=733  石英定压比热 20 degC（与 C^E 差 <1%）
        %   C_E_f=4157 纯水 c_v 20 degC（NIST）
        %   alpha_s=3.45e-5 石英体膨胀 20 degC
        %   alpha_f=2.05e-4 水体膨胀 20 degC（原 7e-6 小 100 倍，已修正）

    case 'chen'      % 陈
        mat.elastic = struct( ...
            'k0',    1e-15,    'morsE', 2,      'visc0', 0.001, ...
            'den_s', 2544,     'poros', 0.203,  'den_f', 940, ...
            'Kf',    2.1e9,    'Ks',    37e9,   ...
            'mu_b',  26.2e9,   'Kb',    20.8e9, ...
            'nu1',   3578e9,   'nu2',  -1310e9, 'nu3',  -345e9, ...
            'gama', -1150e9,   'gama1',  1.93e9, ...
            'gama2',  220e9,   'gama3', -110e9);
        mat.thermal = struct( ...
            'T0',      300,    'k0_s',   12,     'k0_f',   10, ...
            'time0_s', 1e-11,  'time0_f', 1e-9, ...
            'C_E_s',   800,    'C_E_f',  2000, ...
            'K',       1e-7,   'alpha_s', 2e-5,  'alpha_f', 7e-6);

    case 'hou'       % Hou, Fu & Carcione (2024), GEOPHYSICS 89(5) MR297-MR307,
                     %   doi:10.1190/GEO2023-0625.1，Table 1（p.MR299）水饱和列
                     %   注意：与 'zhao' 不是同一块岩石（Kb 7 vs 13.2 GPa、poros 0.1 vs 0.197、
                     %   k0 1 D vs 43.5 mD、morsE 2 vs 2.848；干骨架泊松比 0.189 vs 0.320）
                     %   原文给：Ks=35, Km=7, mu=5.5, Kf=2.4 GPa；den_s=2.65, den_f=1.0 g/cc；
                     %            poros=0.1, kappa=1 D, eta=1e-3 Pa.s, morsE=2；
                     %            热学 T0=400 K, k0_s=12, k0_f=8.76, C_E_s=900, C_E_f=4180,
                     %            alpha_s=2.4e-6, alpha_f=7e-4, time0_s=0.08 ns, time0_f=2.5 ns
                     %   原文未给：三阶常数 nu1..nu3 / gama..gama3  -> 沿用 'zhao'
                     %   原文无相间热交换项（耦合由 F12=F21 非零承担） -> K 沿用原值
                     %   原文另有交叉膨胀 alpha_fs=-2.4e-4，未进程序（本文 Omega 用 alpha_s/alpha_f 显式）
        mat.elastic = struct( ...
            'k0',    1e-12,    'morsE', 2,      'visc0', 2.17e-4, ...
            'den_s', 2650,     'poros', 0.1,    'den_f', 1000, ...
            'Kf',    2.4e9,    'Ks',    35e9,   ...
            'mu_b',  5.5e9,    'Kb',    7e9,    ...
            'nu1',   3578e9,   'nu2',  -1310e9, 'nu3',  -345e9, ...
            'gama', -1150e9,   'gama1',  1.93e9, ...
            'gama2',  220e9,   'gama3', -110e9);
        mat.thermal = struct( ...
            'T0',      400,    'k0_s',   12,     'k0_f',   8.76, ...
            'time0_s', 8e-11,  'time0_f', 2.5e-9, ...
            'C_E_s',   900,    'C_E_f',  4180, ...
            'K',       1e-7,   'alpha_s', 2.4e-6, 'alpha_f', 7e-4);

    case 'zong_sand' % Zong et al. (2023) Surv. Geophys. 44:425-462, Table 2
                     %   「Cylindrical sandstone」那列（水饱和，据 Qu et al. 2018）
                     %   弹性参数全部按原表；热学 Table 2 未给 -> 沿用 'zhao'
                     %   原表：phi=0.004, k0=0.01 md, at=3.7, KM=28.18, Kb=56.46,
                     %         rho_s=2760, mu_b=31.29, Kf=2.25 GPa, rho_f=1000, eta=1e-3
                     %   TOE: nu1'=-39196, nu2'=-9040, nu3'=-1507, gama=-18000,
                     %        gama1=1.93, gama2=21000, gama3=44000  (均 GPa)
        mat.elastic = struct( ...
            'k0',    1e-17,    'morsE', 3.7,    'visc0', 1e-3, ...
            'den_s', 2760,     'poros', 0.004,  'den_f', 1000, ...
            'Kf',    2.25e9,   'Ks',    56.46e9, ...
            'mu_b',  31.29e9,  'Kb',    28.18e9, ...
            'nu1', -39196e9,   'nu2',  -9040e9, 'nu3', -1507e9, ...
            'gama', -18000e9,  'gama1',   1.93e9, ...
            'gama2', 21000e9,  'gama3', 44000e9);
        mat.thermal = struct( ...
            'T0',      300,    'k0_s',   7.69,   'k0_f',   0.598, ...
            'time0_s', 1e-11,  'time0_f', 1e-9, ...
            'C_E_s',   740,    'C_E_f',  4157, ...
            'K',       1e-7,   'alpha_s', 3.45e-5, 'alpha_f', 2.07e-4);

    otherwise
        error('param_0:material', '未知材质预设: %s（可选 zhao / chen / hou / zong_sand）', name);
end
end
