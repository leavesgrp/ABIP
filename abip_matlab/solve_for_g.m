function [g, gTh, itn] = solve_for_g(work, data, h, n, m, k, rho_y, use_indirect, extra_verbose)
    [g, itn]    = solve_lin_sys(work, data, h, n, m, k, zeros(m,1), rho_y, use_indirect, extra_verbose);
    g(m+1:end)  = -g(m+1:end);
    gTh         = g'*h;
end