function run_all()
% run_all -- 按三组 (pc, pf) 批量生成导波图与全波图
%
%   pc=65 MPa, pf=5 MPa
%   pc=45 MPa, pf=5 MPa
%   pc= 0 MPa, pf=0 MPa   （该组无瑞利波）
%
%   输出：重构之后/result/pc<pc>MPa_pf<pf>MPa/ 下
%         fig_<模式>_频散/激发.png|pdf  （n=0 斯通利、n=0 伪瑞利、n=1 弯曲、n=2 螺旋）
%         fig_全波波形_n0/n1/n2.png|pdf

here = fileparts(mfilename('fullpath'));
addpath(here);

cases = [ 65e6,  5e6;
          45e6,  5e6;
           0,    0   ];
% cases = [ 45e6,  5e6];

for ic = 1:size(cases, 1)
    pc = cases(ic,1);  pf = cases(ic,2);
    fprintf('\n========== case %d/%d : pc = %g MPa, pf = %g MPa ==========\n', ...
            ic, size(cases,1), pc/1e6, pf/1e6);
    fprintf('--- guidewave ---\n');
    plot_guidewave(pc, pf);
    % fprintf('--- fullwave ---\n');
    % plot_fullwave(pc, pf, [0 1 2]);
end

fprintf('\nall done.\n');
end