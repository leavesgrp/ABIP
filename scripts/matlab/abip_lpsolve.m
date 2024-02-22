function [x, y, s, info] = abip_lpsolve(data, K, params)
% Invoke ABIP LP solver

A = data.A;
b = data.b;
c = data.c;

% Another check to ensure that no other cones exist
if isfield(K, 'f') || isfield(K, 'q') || isfield(K, 'rq') || ~isfield(K, 'l')
    error("Invalid conic format for LP");
end % End if

lpparams = abipi_lpparam_convert(params);

% Solution method pcg/direct
tic;
if params.pcg
    [x, y, s, lpinfo] = abip_indirect(data, lpparams);
else
    [x, y, s, lpinfo] = abip_direct(data, lpparams);
end % End if
tlp = toc;

info = abipi_lpinfo_convert(lpinfo);
info.time = tlp;

info.pobj = c' * x;
info.dobj = b' * y;

info.solver = 'abip-lp';

end % End function

function [lpparams] = abipi_lpparam_convert(params)
% Convert general interface into ABIP-LP interface
lpparams.verbose = params.verbose;
lpparams.normalize = params.normalize;
lpparams.max_admm_iter = params.max_admm_iter;
lpparams.max_ipm_iters = params.max_ipm_iter;
lpparams.timelimit = params.timelimit;
lpparams.eps = params.tol;

% Convert scaling parameter
lpparams.origin_rescale = 0;
lpparams.pc_ruiz_rescale = 0;
lpparams.qp_rescale = 0;
if params.lpalg.scaling_method == 1
    lpparams.pc_ruiz_rescale = 1;
elseif params.lpalg.scaling_method == 2
    lpparams.qp_rescale = 1;
elseif params.lpalg.scaling_method == 3
    lpparams.origin_rescale = 1;
end % End if

% Convert the rest of parameters
lpparams.restart_thresh = params.lpalg.restart_thresh;
lpparams.restart_freq = params.lpalg.restart_freq;
lpparams.feasopt = params.lpalg.feasopt;
lpparams.half_update = params.lpalg.half_update;

end % End function

function [info] = abipi_lpinfo_convert(lpinfo)
% Convert LP solution statistics to general interface
info.status = lpinfo.status;
info.ipm_iter = lpinfo.ipm_iter;
info.admm_iter = lpinfo.admm_iter;
info.pres = lpinfo.resPri;
info.dres = lpinfo.resDual;
info.gap = lpinfo.relGap;
end % End function