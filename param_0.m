function param = param_0(pc, pf, dT, n, r)
% param_0 —— 不随频率变化的全部参数（材质 + 派生量 + 随 pc 变的静态量 + 几何）
%
% 五个算例量都可以省略，省略时用下面这组默认值（main 里直接 param_0() 即可）：
%   pc      围压 [Pa]        默认 65e6
%   pf      孔压 [Pa]        默认 5e6
%   dT      温度变化 [K]     默认 0
%   n       方位序数（Bessel 阶） 默认 0
%   r       接收器径向位置 [m]    默认 = R（井壁）；R 固定 0.1 m，不作输入
if nargin < 1 || isempty(pc),     pc     = 65e6;    end
if nargin < 2 || isempty(pf),     pf     = 5e6;    end
if nargin < 3 || isempty(dT),     dT     = 0;    end
if nargin < 4 || isempty(n),      n      = 0;    end
R = 0.1;                                        % 井半径 [m]（固定，仅内部用）
if nargin < 5 || isempty(r),      r      = R;   end

% —— 源径（张碧星多极源分布半径；只进多极因子，不影响色散）——
r0 = 0.01;                   % 源分布半径 [m]
s0 = 1e-3;                   % 源强标定（原 fw_src_amp）

% —— 方位角（角向因子 cos(n*(th_rec-th_src))）——
th_rec = 0;                  % 接收方位角 [rad]
th_src = 0;                  % 源方位角 [rad]

%% §4 材质与热学预设（换算例只改这一段）
%% §4.0 缩放因子（材质公式的依赖，故放在本节最前）
scale   = 1;      % 热学分组缩放总倍率：1 = 物理值；1e-1 = 旧口径（热效应整体弱约 10×）
                  %   下游 sc_k / sc_C / sc_t / sc_K / sc_alp 全部取它 —— 改一处即整体缩放热学

%% §4.1 孔隙弹性材质 —— 默认「赵志强」
% —— 渗流（下游进 Biot 耗散项；本文件不直接算它）——
k0    = 43.5e-15;   % 渗透率 [m²]（43.5e-15 m²）
morsE = 2.8481;     % 渗透率动态修正参数（无量纲，沿用原程序口径；本文件派生量不用它）
visc  = 0.00106;    % 流体黏度 [Pa·s]
% —— 密度与孔隙度 ——
den_s = 2581.57;    % 固体颗粒密度 [kg/m³]
poros = 0.197;   % 孔隙度 [-]（基准 0.197 × chengzi）
den_f = 1000;       % 流体密度 [kg/m³]
% —— 模量（§6 的 Alpha / M / lambda_c 全从这里推）——
Kf   = 2.25e9;      % 流体体积模量 [Pa]
Ks   = 37e9;        % 固体颗粒体积模量 [Pa]
mu_b = 5.4e9;       % 干岩石剪切模量 [Pa]
Kb   = 13.2e9;      % 干岩石体积模量 [Pa]

% ---- 备选预设（整段换成对应几行即可）----
% %陈
% k0 = 1e-15; morsE = 2; visc = 0.001;
% den_s = 2544; poros = 0.203; den_f = 940;
% Kf = 2.1e9; Ks = 37e9; mu_b = 26.2e9; Kb = 20.8e9;
%
% %第三套
% k0 = 10e-15; morsE = 3.7; visc = 0.001;
% den_s = 2760; poros = 0.12; den_f = 1000;
% Kf = 2.25e9; Ks = 20.46e9; mu_b = 31.29e9; Kb = 46.18e9;

%% §4.2 VTI 各向异性 —— 默认「赵」
% 两组三阶弹性常数 [Pa]，命名沿用原程序。
% 用途：§8 围压修正用它们把形变 (e_s, zeta_s) 换成 H'/C'/G'/Upsilon'。
nu1 = 3578e9;  nu2 = -1310e9;  nu3 = -345e9;                    % 一组三阶常数
gama = -1150e9;  gama1 = 1.93e9;  gama2 = 220e9;  gama3 = -110e9;  % 另一组三阶常数

% %各向同性（全零）
% nu1 = 0e9; nu2 = 0e9; nu3 = 0e9;
% gama = 0e9; gama1 = 0e9; gama2 = 0e9; gama3 = 0e9;

%% §4.3 热学与耗散（分组缩放：各组倍率可分别调）
% 五个 sc_* 默认都等于 scale；要单独调某一类热学量，就改对应的 sc_*
%   sc_k  ：热导率  ×  → k0_s / k0_f          （越大导热越快）
%   sc_C  ：比热容  ÷  → C_E_s / C_E_f        （注意是"除"：越大热容越小）
%   sc_t  ：弛豫时间 ×  → time0_s / time0_f    （热惯性时间尺度）
%   sc_K  ：孔隙导热 ×  → K                    （热-孔耦合强度）
%   sc_alp：热膨胀  ×  → alpha_s / alpha_f    （热应变驱动强度）
sc_k = scale;  sc_C = scale;  sc_t = scale;  sc_K = scale;  sc_alp = scale;

% —— 温度 ——
T0      = 300;              % 参考温度 [K]（§6 里 F_pij = F_ij/T0 的分母）
                            %   0 = 无温度变化（热膨胀项被乘掉）；改这一行 → theta_s_s/theta_f_s 同步
% —— 导热（× sc_k）——
k0_s    = 1.2e1  * sc_k;    % 固体热导率 [W/(m·K)]
k0_f    = 1e1    * sc_k;    % 流体热导率 [W/(m·K)]
% —— 热弛豫时间（× sc_t）——
time0_s = 1e-11  * sc_t;    % 固体热弛豫时间 [s]
time0_f = 1e-9   * sc_t;    % 流体热弛豫时间 [s]
% —— 比热容（÷ sc_C）——
C_E_s   = 800    / sc_C;    % 固体比热容 [J/(kg·K)]
C_E_f   = 2000   / sc_C;    % 流体比热容 [J/(kg·K)]
% —— 热耦合与热膨胀 ——
K       = 1e-7   * sc_K;    % 热耦合系数（孔隙导热项，沿用原程序量纲）
alpha_s = 2e-5   * sc_alp;  % 固体热膨胀系数 [1/K]
alpha_f = 7e-6   * sc_alp;  % 流体热膨胀系数 [1/K]

%% §5 井外介质派生量（几何参数 R 由主程序给）
% —— Biot 孔弹性模量 ——
Alpha    = 1 - Kb/Ks;                          % Biot 有效应力系数 α（0 < α < 1）
M        = 1/((Alpha-poros)/Ks + poros/Kf);    % Biot 模量 M [Pa]（孔压-体积耦合）
lambda_b = Kb - 2/3*mu_b;                      % 干岩石 Lamé 第一参数 λ_b [Pa]
lambda_c = lambda_b + M*Alpha^2;               % 不排水/封闭 λ_c [Pa]（进 A3333、G_p、C_p）
Kc       = Kb + M*Alpha^2;                     % 封闭体积模量 Kc [Pa]（只用于 Ω_1s / Ω_2s）

% —— 相密度 ——
den_s_prime = den_s*(1 - poros);   % 固相等效密度 (1-φ)ρs [kg/m³]
den_f_prime = den_f*poros;         % 流相等效密度 φρf [kg/m³]
den = 2270;      % 总密度：原程序口径固定 2270（den_s_prime + den_f_prime 算出来也不用）

% —— 热惯性系数（= 等效密度 × 比热容）——
F11 = den_s_prime*C_E_s;   % 固相热惯性；同时等于 -2*G*T0
F22 = den_f_prime*C_E_f;   % 流相热惯性；同时等于 -2*H*T0
F12 = 0;  F21 = 0;         % 本项目口径：两相热惯性不交叉

% —— 热-孔弹性耦合系数（温度变化 → 应力 / 孔压）——
Omega_1s = poros*M*(1-poros)/(1-Alpha)*(Kc/poros/M-Alpha)*alpha_s;   % 固相温度 → 总应力
Omega_1f = poros*M*alpha_f*Alpha;                                    % 流相温度 → 总应力
Omega_2s = poros*M*(1-poros)/(1-Alpha)*(Alpha/poros-1)*alpha_s;      % 固相温度 → 孔压
Omega_2f = poros*M*alpha_f;                                          % 流相温度 → 孔压

% —— 温度变化与热惯性（除以 T0，直接进 8×8 矩阵）——
theta_s_s = dT;  theta_f_s = dT;   % 固/流两相取同一温度变化 ΔT
F_p11 = F11/T0;  F_p12 = F12/T0;  F_p21 = F21/T0;  F_p22 = F22/T0;

% —— 相热导率 ——
k_f = poros*k0_f;  k_s = (1-poros)*k0_s;   % 流相 / 固相热导率 [W/(m·K)]

% —— 井内流体（井内 2 支模式用；口径：Knf = Kf，不引入 Q 衰减）——
Knf = Kf;  k_f_in = k0_f;                  % 井内流体体积模量 / 热导率
Omega_2f_in = alpha_f*Knf;                 % 井内流体热膨胀耦合
den_f_in = den_f;                          % 井内流体密度
F_p22_in = den_f_in*C_E_f/T0;              % 井内流体热惯性（÷T0）
time0_f_in = time0_f;                      % 井内流体热弛豫时间



% —— 驱动量（有效围压 / 有效孔压，含热应变项）——
Qiues    = -pc + Omega_1s*theta_s_s + Omega_1f*theta_f_s;   % 作用在固相上的有效围压
Qiuzetas =  pf - Omega_2s*theta_s_s - Omega_2f*theta_f_s;   % 作用在孔隙上的有效孔压

% —— 静态形变 ——
e_s    = (Qiues + Alpha*Qiuzetas)/Kb;                       % 骨架体积应变
zeta_s = (Alpha*Qiues + (lambda_c/M)*Qiuzetas)/Kb;          % 孔隙体积应变增量
e33    = e_s/3;                                             % = 原程序的 u_s_33（耦合项用它）
w33    = zeta_s/3;                                          % = 原程序的 w_s_33
AlphaM = Alpha*M;                                           % αM 在 A3333/G'/C'/Υ 里反复出现

% —— 四个模量（干岩石模量 + 三阶常数修正）——
A3333 = lambda_c + 2*mu_b ...
        + (lambda_c+nu1+2*nu2)*e_s ...
        + 2*(lambda_c+3*mu_b+2*nu2+4*nu3)*e33 ...
        - zeta_s*(AlphaM-2*gama2);                          % 轴向模量 C3333 [Pa]
G_p = mu_b + (lambda_c + nu2 + (4/3)*(mu_b + nu3))*e_s ...
        - (AlphaM + gama/2)*zeta_s;                         % 修正剪切模量 G' [Pa]
C_p = AlphaM + (AlphaM+gama)*e33 ...
        - (2*gama2+gama)*e_s ...
        - zeta_s*(AlphaM+2*gama3);                          % 修正耦合模量 C' [Pa]
Upsilon = M + e_s*(2*gama3+AlphaM) ...
        - zeta_s*3*(M-2*gama1);                             % 修正孔压模量 Υ [Pa]

% —— 再扣掉热应力项，得到最终 H' / Upsilon' ——
H_p       = A3333   - (Omega_1f*theta_f_s + Omega_1s*theta_s_s);
Upsilon_p = Upsilon - (Omega_2f*theta_f_s + Omega_2s*theta_s_s);

% —— 动态耦合系数（带 (1+形变) 因子）——
Omega_11 = Omega_1s*(1 + e33);  Omega_12 = Omega_1f*(1 + e33);
Omega_21 = Omega_2s*(1 + w33);  Omega_22 = Omega_2f*(1 + w33);


% —— 打包（字段顺序 = §8 并入 param 的顺序，勿随意调整）——
param = struct( ...
    'A3333',A3333, 'G_p',G_p, 'C_p',C_p, ...
    'Upsilon',Upsilon, 'H_p',H_p, 'Upsilon_p',Upsilon_p, ...
    'Omega_11',Omega_11, 'Omega_12',Omega_12, ...
    'Omega_21',Omega_21, 'Omega_22',Omega_22, ...
    'F_p11',F_p11, 'F_p12',F_p12, 'F_p21',F_p21, 'F_p22',F_p22, ...
    'k_f',k_f, 'k_s',k_s, ...
    'time0_s',time0_s, 'time0_f',time0_f, 'T0',T0, 'K',K, ...
    'den',den, 'den_f',den_f, 'poros',poros, ...
    'k0',k0, 'morsE',morsE, 'visc',visc, ...
    'k_f_in',k_f_in, 'Knf',Knf, 'Omega_2f_in',Omega_2f_in, ...
    'den_f_in',den_f_in, 'F_p22_in',F_p22_in, 'time0_f_in',time0_f_in, ...
    'R',R, 'n',n, 'r0',r0, 's0',s0, 'r',r, 'dT',dT, 'th_rec',th_rec, 'th_src',th_src, 'pc',pc);
end
