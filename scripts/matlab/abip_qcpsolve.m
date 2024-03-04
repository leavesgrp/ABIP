function [x, y, s, info] = abip_qcpsolve(data, K, params)
% Invoke ABIP QCP solver

% Check if the data format is more suitable for LP
if ~isfield(K, 'f') && ~isfield(K, 'q') && ~isfield(K, 'rq') && ~isfield(K, 'z')
    warning("The problem is an LP. It is recommended to use the LP solver");
end % End if

qcpparams = abipi_qcpparam_convert(params);

tic;
[sol, qcpinfo] = abip_qcp(data, K, qcpparams);
tqcp = toc;

x = sol.x;
y = sol.y;
s = sol.s;

info = abipi_qcpinfo_convert(qcpinfo);
info.time = tqcp;

info.solver = 'abip-qcp';

end % End function

function [qcpparams] = abipi_qcpparam_convert(params)
% Convert general interface into ABIP-QCP interface
qcpparams.verbose = params.verbose;
qcpparams.normalize = params.normalize;
qcpparams.max_admm_iter = params.max_admm_iter;
qcpparams.max_ipm_iters = params.max_ipm_iter;
qcpparams.timelimit = params.timelimit;

% Convert linear system
if params.pcg
    qcpparams.linsys_solver = 3;
else
    % There are other options for indirect solver, but we omit them in this
    % version
    qcpparams.linsys_solver = 1;
end % End if

% Convert tolerances
qcpparams.eps_p = params.tol;
qcpparams.eps_d = params.tol;
qcpparams.eps_g = params.tol;

% Convert the rest
qcpparams.rho_x = params.qcpalg.rho_primal;
qcpparams.rho_y = params.qcpalg.rho_dual;
qcpparams.psi = params.qcpalg.admm_tol_factor;

end % End function

function [info] = abipi_qcpinfo_convert(qcpinfo)
% Convert QCP solution statistics to general interface
info.status = qcpinfo.status;
info.ipm_iter = qcpinfo.ipm_iter;
info.admm_iter = qcpinfo.admm_iter;
info.pres = qcpinfo.res_pri;
info.dres = qcpinfo.res_dual;
info.gap = qcpinfo.gap;
info.pobj = qcpinfo.pobj;
info.dobj = qcpinfo.dobj;
info.gap = qcpinfo.gap;
end % End function