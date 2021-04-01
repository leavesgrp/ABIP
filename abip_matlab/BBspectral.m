function beta = BBspectral(work, data, mu, n, m, l, k, u, v, rho_y, use_indirect, extra_verbose, ...
    h, g, gTh, beta, alpha, eps_cor, eps_pen)
    %% intialization
    beta_pre = beta;
    
    %% the first iterate
    u_pre = u;  
    v_pre = v; 
    
    
    %% main loop
    for j=1:5,
        %% the second iterate
        [ut, ~]     = project_lin_sys(work, data, n, m, k, u_pre, v_pre, rho_y, use_indirect, ...
            extra_verbose, h, g, gTh);
        rel_ut      = alpha*ut+(1-alpha)*u_pre;
        rel_ut(1:m) = ut(1:m);                       
        u           = rel_ut - v_pre;
        temp        = u(m+1:end);
        u(m+1:end)  = 0.5*(temp+sqrt(temp.*temp+4*mu/beta_pre));
        v           = v_pre + (u - rel_ut);
        
        %% the third iterate
        [ut_next, ~]     = project_lin_sys(work, data, n, m, k, u, v, rho_y, use_indirect, ...
            extra_verbose, h, g, gTh);
        rel_ut_next      = alpha*ut_next+(1-alpha)*u;
        rel_ut_next(1:m) = ut_next(1:m);                       
        u_next           = rel_ut_next - v;
        temp             = u_next(m+1:end);
        u_next(m+1:end)  = 0.5*(temp+sqrt(temp.*temp+4*mu/beta_pre));
        v_next           = v + (u_next - rel_ut_next);
        
        %% curvature information
        v_gap    = v_next-v-(1-alpha)*(u_next-u);
        ut_gap   = v+u_next-v_next-v_pre-u+v; 
        u_gap    = u-u_next; 
        vv       = v_gap'*v_gap;
        vut      = v_gap'*ut_gap;
        utut     = ut_gap'*ut_gap;
        vu       = v_gap'*u_gap;
        uu       = u_gap'*u_gap;
        v_norm   = norm(v_gap);
        ut_norm  = norm(ut_gap);
        u_norm   = norm(u_gap);
        
        %% SD and MG spectral stepsize
        alpha_SD = vv/vut; 
        alpha_MG = vut/utut;
        gamma_SD = vv/vu; 
        gamma_MG = vu/uu;
        % fprintf('alpha_SD: %3.6f, alpha_MG: %3.6f, gamma_SD: %3.6f, gamma_MD: %3.6f\n', alpha_SD, alpha_MG, gamma_SD, gamma_MG);
        if 2*alpha_MG > alpha_SD
            alpha_ss = alpha_MG; 
        else
            alpha_ss = alpha_SD - 0.5*alpha_MG; 
        end
        
        if 2*gamma_MG > gamma_SD
            gamma_ss = gamma_MG; 
        else
            gamma_ss = gamma_SD - 0.5*gamma_MG; 
        end
        
        %% safeguarding. 
        alpha_cor = vut/(v_norm*ut_norm);
        gamma_cor = vu/(v_norm*u_norm); 
        
        %% spectral stepsize
        if alpha_cor>eps_cor && gamma_cor>eps_cor
            beta = sqrt(alpha_ss*gamma_ss); 
        elseif alpha_cor>eps_cor && gamma_cor<=eps_cor 
            beta = alpha_ss; 
        elseif alpha_cor<=eps_cor && gamma_cor>eps_cor 
            beta = gamma_ss; 
        else 
            beta = beta_pre; 
        end
        
        if abs(beta-beta_pre)>0 && abs(beta-beta_pre)<=eps_pen
            beta = 0.5*(beta+beta_pre); 
            break; 
        elseif abs(beta-beta_pre)>eps_pen
            beta_pre = beta;
            u_pre = u; 
            v_pre = v; 
            v_pre(m+1:l) = (mu/beta_pre)./u_pre(m+1:l);
        else
            u_pre = u; 
            v_pre = v; 
        end
        
        % fprintf('beta: %3.7e, vv: %3.7e, uu: %3.7e, vut: %3.7e, vu: %3.7e, utut: %3.7e\n', beta, vv, uu, vut, vu, utut);
    end
end