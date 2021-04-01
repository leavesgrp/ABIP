function [x, y, s, info] = abip(varargin)
    % COPL 1.0.0 for version call: abip_version()

    data = varargin{1};
    if nargin >= 2
        pars = varargin{2};
    else
        pars = [];
    end
    
    if (isfield(pars,'use_direct') && pars.use_direct)
    	[x, y, s, info] = abip_direct(data, K, pars);
    elseif (isfield(pars,'gpu') && pars.gpu)
    	[x, y, s, info] = abip_gpu(data, K, pars);
    else
   	[x, y, s, info] = abip_indirect( data, K, pars);
    end