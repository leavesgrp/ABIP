function [y, itn] = solve_lin_sys(w, data, rhs, n, m, k, warm_start, rho_y, use_indirect, extra_verbose)
    % assumes matrix [rho_y*speye(m) data.A; data.A' -speye(n)] 
    if use_indirect
        if k == -1 
            tol = 1e-9 * norm(rhs(1:m)); 
        else
            tol = max(1e-9, 1e-5 * norm(rhs(1:m)) / (k+1)^2);  
        end
        y               = zeros(m+n, 1); 
        [y(1:m), itn]   = pcg_abips(data.A, rhs(1:m)+data.A*rhs(m+1:m+n), warm_start(1:m), w.M, ...
            rho_y, m, tol, extra_verbose);  
        y(m+1:m+n)      = -rhs(m+1:m+n) + data.A'*y(1:m);
    else
        y               = (w.L'\(w.d\(w.L\rhs(w.P))));
        y(w.P)          = y;
        itn             = -1; 
    end
end