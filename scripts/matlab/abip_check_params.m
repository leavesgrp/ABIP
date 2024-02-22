function [params] = abip_check_params(params)
% Check validity of ABIP parameters

default_params = abip_get_params();

% Check common ABIP parameters
if ~isfield(params, 'verbose')
    params.verbose = default_params.verbose;
end % Endif

if ~isfield(params, 'normalize')
    params.normalize = default_params.normalize;
end % Endif

if ~isfield(params, 'pcg')
    params.pcg = default_params.pcg;
end % End if

if ~isfield(params, 'max_admm_iter')
    params.max_admm_iter = default_params.max_admm_iter;
end % End if

if params.max_admm_iter <= 0
    error("Invalid max_admm_iter. Must be > 0");
end % End if

if ~isfield(params, 'max_ipm_iter')
    params.max_ipm_iter = default_params.max_ipm_iter;
end % End if

if params.max_ipm_iter <= 0
    error("Invalid max_ipm_iter. Must be > 0");
end % End if

if ~isfield(params, 'timelimit')
    params.timelimit = default_params.timelimit;
end % End if

if params.timelimit <= 0.0
    error("Invalid timelimit. Must be > 0.0");
end % End if

if ~isfield(params, 'tol')
    params.tol = default_params.tol;
end % End if

if params.tol <= 0.0
    error("Invalid tol. Must be > 0.0");
end % End if

if ~isfield(params, 'lpalg')
    params.lpalg = default_params.lpalg;
end % Endi if
params.lpalg = abipi_check_lpparams(params.lpalg);

if ~isfield(params, 'qcpalg')
    params.qcpalg = default_params.qcpalg;
end % Endi if
params.qcpalg = abipi_check_qcpparams(params.qcpalg);

end % End function

function [params] = abipi_check_lpparams(params)

default_lpparams = abip_get_params();
default_lpparams = default_lpparams.lpalg;

if ~isfield(params, 'restart_thresh')
    params.restart_thresh = default_lpparams.restart_thresh;
end % End if

if params.restart_thresh <= 0
    error("Invalid LP parameter restart_thresh. Must be > 0");
end % End if

if ~isfield(params, 'restart_freq')
    params.restart_freq = default_lpparams.restart_freq;
end % End if

if params.restart_freq <= 0
    error("Invalid LP parameter restart_freq. Must be > 0");
end % End if

if ~isfield(params, 'feasopt')
    params.feasopt = default_lpparams.feasopt;
end % End if

if ~isfield(params, 'scaling_method')
    params.scaling_method = default_lpparams.scaling_method;
end % End if

if ~ismember(params.scaling_method, [1, 2, 3])
    error("Invalid scaling method. Must be one of [1, 2, 3].");
end % End if 

if ~isfield(params, 'half_update')
    params.half_update = default_lpparams.half_update;
end % End if

end % End internal function

function [params] = abipi_check_qcpparams(params)

default_qcpparams = abip_get_params();
default_qcpparams = default_qcpparams.qcpalg;

if ~isfield(params, 'rho_primal')
    params.rho_primal = default_qcpparams.rho_primal;
end % End if

if params.rho_primal <= 0
    error("Invalid QCP parameter rho_primal. Must be > 0");
end % End if

if ~isfield(params, 'rho_dual')
    params.rho_dual = default_qcpparams.rho_dual;
end % End if

if params.rho_dual <= 0
    error("Invalid QCP parameter rho_dual. Must be > 0");
end % End if

if ~isfield(params, 'admm_tol_factor')
    params.admm_tol_factor = default_qcpparams.admm_tol_factor;
end % End if

if params.admm_tol_factor <= 0
    error("Invalid QCP parameter admm_tol_factor. Must be > 0");
end % End if

end % End internal function