function e = excite(fr, param, kz, omega)
% excite —— 极点处的激发留数（张碧星多极源）
%
%   井壁 / 井轴的振幅已移入 M8_solve（第 8、9 输出 an_wall / an_axis），本函数只做叠加：
%       A_wall = A .* I_n(knr*r) .* e^{|Re knr r| - |Re 2knr R|}
%       A_axis = A .* (knr/2)^n .* e^{-|Re 2knr R|}
%       A      = 1i * an .* exp(1i*kz*z)          （z 由 M8_solve 的 z_obs 给）
%
%   输出（两支叠加后的标量）：
%       e.ap_wall = sum(A_wall)              井壁声压激发（待乘 rho_f*omega^2）
%       e.at_wall = sum(A_wall .* alphani)   井壁温度激发
%       e.ap_axis = sum(A_axis)              井轴声压激发
%       e.at_axis = sum(A_axis .* alphani)   井轴温度激发
%       e.ok      留数是否可信（an_o 全有限）


% 一次组装：直接取 an_wall / an_axis（第 8、9 输出）
[~, ~, ~, ~, an_wall, an_axis] = M8_solve(fr, param, kz, 3);   %#ok<ASGLU>

e = struct('ap_wall', [], 'at_wall', [], 'ap_axis', [], 'at_axis', [], 'ok', false);

aln = fr.alphani(:).';                      % 1x2 温度/势比（与 an_wall 同向）

e.ap_wall = sum(an_wall);
e.at_wall = sum(an_wall .* aln);
e.ap_axis = sum(an_axis);
e.at_axis = sum(an_axis .* aln);

e.ok = all(isfinite(an_wall)) && all(isfinite(an_axis));   % 留数是否可信
end