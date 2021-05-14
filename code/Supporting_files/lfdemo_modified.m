% ldfwdemo.m
% m file to demo loadflow solutions
% a loadflow report is placed in lf_rep.txt if requested
% clearvars -except dfile
% clear global
% dfile = 'data16mt9eg18_2.m';
% global variables
pst_var 
global gen_chg_idx

disp('loadflow demo program')

lfile =length(dfile);
  
dfile = lower(dfile(1:lfile-2));
eval(dfile);

nbus = length(bus(:,1));
gen_chg_idx = ones(nbus,1);
[bus_sol,line_sol,line_flow] = loadflow(bus,line,1e-9,10,1.0,'y',2);
 
