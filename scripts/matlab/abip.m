function [x, y, s, info] = abip(data, K, params)
% ABIP: An ADMM-based Interior Point Method
% TODO: More details here

if nargin == 0
    fprintf("ABIP: An ADMM-based Interior Point Method \n");
    fprintf("Usage: [x, y, s] = abip(data, K, params \n");
end % End if

if nargin < 2
    error("Invalid number of inputs. Expected at least data and K.");
end % End if

if nargin == 2
    params = abip_get_params();
end % End if 

% Do some parameter check
params = abip_check_params(params);

% TODO: Add conversion for linear constraints if 'f' exists
if isfield(K, 'f') || isfield(K, 'q') || isfield(K, 'rq') || params.solver == 1
    % Call QCP Solver
    [x, y, s, info] = abip_qcpsolve(data, K, params);
else
    % Call LP Solver 
    [x, y, s, info] = abip_lpsolve(data, K, params);
end % End if

end % End function