function [ut, itn] = project_lin_sys(work, data, n, m, k, u, v, rho_y, use_indirect, extra_verbose, h, g, gTh)    
    ut                  = u+v;
    ut(1:m)             = rho_y*ut(1:m);
    ut(1:m+n)           = ut(1:m+n) - ut(end)*h;
    ut(1:m+n)           = ut(1:m+n) - h*((g'*ut(1:m+n))/(gTh+1));
    warm_start          = u(1:n+m);
    ut(m+1:end-1)       = -ut(m+1:end-1);
    [ut(1:m+n), itn]    = solve_lin_sys(work, data, ut(1:m+n), n, m, k, warm_start, rho_y, ...
        use_indirect, extra_verbose);
    ut(end)             = (ut(end) + h'*ut(1:m+n));    
end