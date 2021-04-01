function [x, objp] = postsolve(x, info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postsolve recover the original variables for the solution.                     
% .
% Usage:    [x, objp] = postsolve(x, info)
% Problem:   min c'x, s.t. Ax=b, lbounds<=x<=ubounds. 
% Input:    Variable:        x.
%           Information:     info.
% 
% Author: Tianyi Lin, UC Berkeley
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% recover the original data
c = info.c; 
lbounds = info.lbounds;

%% recover the statue
fixed_exist  = info.fixed_exist; 
zrcols_exist = info.zrcols_exist;
sgtons_exist = info.sgtons_exist;
free_exist   = info.free_exist;
lbounds_non0 = info.lbounds_non0;
n            = info.n; 

%% recover the reformulated LP.
x = x(1:n); 

%% recover the lower bound.
if (lbounds_non0)
    x = x + lbounds; 
end;

%% recover the free variables. 
if (free_exist)
    ifree         = info.ifree; 
    infree        = info.infree; 
    nfree         = info.nfree;
    x_minus       = x(n-nfree+1:n); 
    x_plus        = x(n-2*nfree+1:n-nfree);
    tmp(infree)   = x(1:n-2*nfree);
    tmp(ifree)    = x_plus - x_minus; 
    x             = sparse(tmp);      
end

%% recover singleton rows.
if (sgtons_exist)
    isolved       = info.isolved; 
    insolved      = info.insolved; 
    xsolved       = info.xsolved;
    tmp(insolved) = x;
    tmp(isolved)  = xsolved;
    x             = sparse(tmp);
end;

%% recover zero columns.
if (zrcols_exist)
    izrcol      = info.izrcol;  
    inzcol      = info.inzcol;
    xzrcol      = info.xzrcol;
    tmp(inzcol) = x;
    tmp(izrcol) = xzrcol;
    x           = sparse(tmp);
end;


%% recover fixed variables.
if (fixed_exist)
    ifix      = info.ifix; 
    infx      = info.infx; 
    xfix      = info.xfix;
    tmp(infx) = x; 
    tmp(ifix) = xfix;
    x         = sparse(tmp);
end;

if (size(x,1) < size(x,2))
   x = x';
end;

objp = full(c'*x);