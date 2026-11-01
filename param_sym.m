function [C_sym, jie_sym, det_sym, vars] = param_sym()
% param_sym —— 只做符号推导，输出三个 sym 表达式（结果 persistent 缓存，只算一次）
%
%   M0(S_sq) = [ S_sq*A_p - den_m,      O    ;
%                S_sq*O',          S_sq*K - F ]
%
%   C_sym   = det(M0(1:4,1:4)) 关于 S_sq 的系数 [C4 C3 C2 C1 C0]
%   jie_sym = -adj(M0(2:4,2:4)) * M0(2:4,1)      （alpha 的分子）
%   det_sym = det(M0(2:4,2:4))                    （alpha 的分母）
%   vars    = {S_sq, H_p, C_p, Upsilon_p, Ks~, Kf~, O11, O12, O21, O22,
%              den, den_f, den_t, F11, F12, F21, F22}
%
% get_poly_coeffs 用 vars(2:end)；get_jie/get_det 用 vars。

persistent C J D V
if ~isempty(C)
    C_sym = C;  jie_sym = J;  det_sym = D;  vars = V;
    return
end

syms S_sq sym_Ks sym_Kf sym_den_t
syms sym_H_p sym_C_p sym_Upsilon_p sym_O11 sym_O12 sym_O21 sym_O22 real
syms sym_den sym_den_f sym_F11 sym_F12 sym_F21 sym_F22 real

A_p   = [sym_H_p, sym_C_p;  sym_C_p, sym_Upsilon_p];
den_m = [sym_den, sym_den_f; sym_den_f, sym_den_t];
O     = [sym_O11, sym_O12;  sym_O21, sym_O22];
F     = [sym_F11, sym_F12;  sym_F21, sym_F22];
K     = [sym_Ks, 0;         0, sym_Kf];

M0 = [S_sq*A_p - den_m, O;  S_sq*O', S_sq*K - F];

C_sym   = coeffs(det(M0(1:4,1:4)), S_sq, 'All');
jie_sym = -adjoint(M0(2:4,2:4)) * M0(2:4,1);
det_sym = det(M0(2:4,2:4));
vars    = {S_sq, sym_H_p, sym_C_p, sym_Upsilon_p, sym_Ks, sym_Kf, ...
           sym_O11, sym_O12, sym_O21, sym_O22, ...
           sym_den, sym_den_f, sym_den_t, sym_F11, sym_F12, sym_F21, sym_F22};

C = C_sym;  J = jie_sym;  D = det_sym;  V = vars;
end
