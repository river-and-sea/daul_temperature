function [kz, kzpath] = kz_newton(detf, kz0, omega)
% kz_newton —— 割线牛顿（原程序口径）
%   detf   @(kz) 复数行列式；内部对 |det| 求极小
%   kz0    初值（粗扫值或外推值）
%
% 两个方向的复步长交替、步长衰减 jianban、缩放 suofangbi，逐句照抄原程序。
% 返回 kzpath：每轮迭代"更新前"的 kz（最后一轮即接受点）。

detabs = @(z) abs(detf(z));

cdzk0   = 1e-6;           % 割线复步长系数（方向 1）
cdzk1   = 1i*1e-6;        % 割线复步长系数（方向 2，与方向 1 差 90°）
amin    = 1e-16;          % 收敛阈值系数：|dkz| <= omega*amin 才停
maxu    = 64;             % 最大迭代步数
step_bi = 0.95;           % 步长衰减因子

dk1 = cdzk0 * kz0^2/omega;
dk0 = cdzk1 * kz0^2/omega;

kz        = kz0;
kzpath    = nan(1, maxu);
kk        = 0;
jianban   = 1;
dkz       = 1;
suofangbi = 1;
DD_pre    = Inf;
jio       = 1;

while abs(dkz) > omega*amin
    kk = kk + 1;

    if jio
        cdzk = dk1;
    else
        cdzk = dk0;
    end
    cdzk = cdzk*jianban;

    DD1 = detabs(kz);
    if DD_pre < DD1
        DD1    = DD_pre;
        kz     = kzpath(kk-1);
        jianban = jianban*0.8;
        jio    = ~jio;
    else
        DD_pre = DD1;
    end

    DD2     = detabs(kz + cdzk);
    dDdk    = (DD2 - DD1) / cdzk;      % |det| 的割线斜率（最后一步即 D' 的估计）
    dkz     = -DD1/dDdk*suofangbi;
    kzpath(kk) = kz;

    suofangbi = suofangbi*step_bi;
    kz        = kz + dkz;

    if kk >= maxu
        break
    end
end

kz = kzpath(kk);          % 最后一轮"更新前"的值就是接受点
end
