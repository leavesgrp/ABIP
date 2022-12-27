function [params] = abip_get_params()
% Set the default parameters for ABIP: An ADMM-based Interior Point Method
% TODO: Write a brief document here.

params.verbose       = 1;
params.normalize     = 1;
params.pcg           = 0;
params.max_admm_iter = 1000000;
params.max_ipm_iter  = 500;
params.timelimit     = 3600;
params.tol           = 1e-03;
params.solver        = -1;

params.lpalg  = abipi_get_lp_params();
params.qcpalg = abipi_get_qcp_params();

end % End function

function [lpalg] = abipi_get_lp_params()

lpalg.restart_thresh = 100000;
lpalg.restart_freq   = 1000;
lpalg.feasopt        = 0;
lpalg.scaling_method = 1;
lpalg.half_update    = 0;

end % End internal function

function [qcpalg] = abipi_get_qcp_params()

qcpalg.rho_primal      = 1.0;
qcpalg.rho_dual        = 1e-06;
qcpalg.admm_tol_factor = 1.0;

end % End internal function