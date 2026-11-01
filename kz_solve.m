function [kz, kz_pre] = kz_solve(detf, omega, omega_pre, v, kz_pre)
% kz_solve —— 一个频点上一个井孔模式：首点粗扫 / 其余外推，然后牛顿
%   detf     @(kz) 复数行列式（fr 在闭包里）
%   v        粗扫相速度区间 [v1, v2]（v2 通常取横波波数，原程序口径）
%   kz_pre   上一频点的 kz（首点传 NaN）

if isempty(kz_pre) || ~isfinite(kz_pre)
    % 首点：粗扫，直接取 |det| 最小的一点作初值
    att_lo = 1e-6;        % 衰减窗 Im/Re 下限（对数分布，须 > 0）
    att_hi = 1e-1;        % 衰减窗 Im/Re 上限
    Nv     = 200;         % 相速度方向点数
    Na     = 40;          % 衰减方向点数

    vg  = linspace(v(1), v(2), Nv);
    att = logspace(log10(att_lo), log10(att_hi), Na);
    KZ  = (omega ./ vg) .* (1 + 1i*att(:));       % Na x Nv，虚部 = 衰减 * 实部
    F   = arrayfun(@(z) abs(detf(z)), KZ);
    [~, i] = min(F(:));
    kz0 = KZ(i);
else
    kz0 = kz_pre * (omega / omega_pre);          % 原程序的外推规则
end

kz     = kz_newton(detf, kz0, omega);
kz_pre = kz;                                      % 下一频点从这里外推
end
