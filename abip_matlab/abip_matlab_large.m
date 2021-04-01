function [x,y,s,info] = abip_matlab_large(data,params)
%------------------------------------------------------------------------
% linear programming solver, solves:
%
% min c'x
% subject to Ax = b, x >= 0. 
%
% where x \in R^n. 
%
% data: A, b, c, where A,b,c used as above
% 
% params: 
% max_outiter: maximum num iterations for ipm. 
% max_iters:   maximum num iterations for admm.
% eps:         quitting tolerances. 
% alpha:       relaxation parameter (alpha = 1 is unrelaxed).
% sigma:       aggressiveness measure from the path-following scheme. 
% normalize:   heuristic normalization procedure. 
% scale:       heuristic re-scaline procedure.
% rho_y:       y equality rescaling. 
% adaptive:    heuristic barzilai-borwein spectral procedure. 
% eps_cor:     curvature information safeguarding. 
%
% Author: Tianyi Lin, UC Berkeley, 2018.  
%--------------------------------------------------------------------------

%% default settings
max_outiters = 100; 
max_iters    = 1000000;             
eps          = 1e-3;                   
alpha        = 1.8; 
mu           = 1.0; 
normalize    = 1;                
scale        = 1;  
rho_y        = 1e-3;
sigma        = 0.3; 
adaptive     = 1;                 
eps_cor      = 0.2;
eps_pen      = 0.1; 

% conjugate gradient (CG) settings:
use_indirect    = true;   % use conjugate gradient rather than direct method
extra_verbose   = false;  % CG prints summary

%% constants
undet_tol    = 1e-18;             % tol for undetermined solution (tau = kappa = 0)

%% parameter setting 
if nargin==2
    if isfield(params,'max_ipm_iters');   max_outiters = params.max_ipm_iters;   end
    if isfield(params,'max_admm_iters');  max_iters    = params.max_admm_iters;      end
    if isfield(params,'eps');             eps          = params.eps;            end
    if isfield(params,'alpha');           alpha        = params.alpha;          end
    if isfield(params,'sigma');           sigma        = params.sigma;          end
    if isfield(params,'normalize');       normalize    = params.normalize;      end
    if isfield(params,'scale');           scale        = params.scale;          end
    if isfield(params,'rho_y');           rho_y        = params.rho_y;          end
    if isfield(params,'adaptive');        adaptive     = params.adaptive;       end
    if isfield(params,'eps_cor');         eps_cor      = params.eps_cor;        end  
    if isfield(params,'eps_pen');         eps_pen      = params.eps_pen;        end   
end

%% data setting 

% dimension
n = length(data.c); 
m = length(data.b); 
l = n+m+1;
u = zeros(l, 1);
v = zeros(l, 1);

nm_c = norm(data.c); 
nm_b = norm(data.b);

%% data normalization
work = struct();
if (normalize)
    [data, work] = normalize_data(data, scale, work);
    D = work.D;
    E = work.E;
    sc_b = work.sc_b;
    sc_c = work.sc_c;
else
    scale = 1;
    D = ones(m,1);
    E = ones(n,1);
    sc_b = 1;
    sc_c = 1;
end

% the matrix Q as in paper
Q = sparse([zeros(m) data.A -data.b; -data.A' zeros(n) data.c; data.b' -data.c' 0]);

%%
if use_indirect
    work.M = 1 ./ diag(rho_y*speye(m) + data.A*data.A'); % pre-conditioner
else
    W                           = sparse([rho_y*speye(m) data.A; data.A' -speye(n)]);  
    [work.L, work.d, work.P]    = ldl(W, 'vector');
end

h           = [-data.b; data.c];
[g, gTh, ~] = solve_for_g(work, data, h, n, m, -1, rho_y, use_indirect, extra_verbose);

gamma = 2;
final_check = 0;
% double_check = 0;

beta = 1;
err_pri = 0; 
err_dual = 0;
gap = 0;

%% Initialization
u(m+1:l) = ones(n+1,1)*sqrt(mu/beta); 
v(m+1:l) = ones(n+1,1)*sqrt(mu/beta);
k        = 0; 

tic;
for i=0:max_outiters-1   
    for j=0:1/(mu)^(1/3)
        % u_pre = u; 
        % v_pre = v;
        
        %% solve linear system
        [ut, ~] = project_lin_sys(work, data, n, m, k, u, v, rho_y, use_indirect, extra_verbose, ...
            h, g, gTh);
    
        %% solve barrier subproblem
        rel_ut      = alpha*ut+(1-alpha)*u;
        rel_ut(1:m) = ut(1:m);                       
        u           = rel_ut - v;
        temp        = u(m+1:end)/2;
        u(m+1:end)  = temp+sqrt(temp.*temp+mu/beta);
    
        %% dual update:
        v = v + (u - rel_ut);
    
        %% convergence checking:
        % err_inner = norm([u-u_pre; v-v_pre])/(1+norm([u;v])+norm([u_pre;v_pre]));
        err_inner = norm(Q*u-v)/(1+norm([u;v]));
        tol = gamma*mu;
        k = k+1;
        % fprintf('||Qu-v|| = %3.6f, gamma = %3.6f, mu = %3.6f\n', err_inner, gamma, mu);
        if err_inner < tol
            break;
        end
        if (final_check && mod(j+1,1)==0)
            tau = abs(u(end));
            kap = abs(v(end)) / (sc_b * sc_c * scale);
            y   = u(1:m) / tau;
            x   = u(m+1:m+n) / tau;
            s   = v(m+1:m+n) / tau;
    
            err_pri  = norm(D.*(data.A * x - data.b)) / (1 + nm_b) / (sc_b * scale); 
            err_dual = norm(E.*(data.A' * y + s - data.c)) / (1 + nm_c) / (sc_c * scale); 
            pobj     = data.c' * x / (sc_c * sc_b * scale);
            dobj     = data.b' * y / (sc_c * sc_b * scale);
            gap      = abs(pobj - dobj) / (1 + abs(pobj) + abs(dobj));
            
            error_ratio = max(gap,max(err_pri,err_dual))/eps;
            solved = error_ratio < 1;
        
            if solved; 
                break; 
            end
            
            if k+1>max_iters
                break; 
            end
        end
    end
    
    %% convergence checking:
    tau = abs(u(end));
    kap = abs(v(end)) / (sc_b * sc_c * scale);
    y   = u(1:m) / tau;
    x   = u(m+1:m+n) / tau;
    s   = v(m+1:m+n) / tau;
    err_pri  = norm(D.*(data.A * x - data.b)) / (1 + nm_b) / (sc_b * scale); 
    err_dual = norm(E.*(data.A' * y + s - data.c)) / (1 + nm_c) / (sc_c * scale); 
    pobj     = data.c' * x / (sc_c * sc_b * scale);
    dobj     = data.b' * y / (sc_c * sc_b * scale);
    gap      = abs(pobj - dobj) / (1 + abs(pobj) + abs(dobj));
    
    if (data.c'*u(m+1:m+n) < 0)
        unb_res = norm(E.*data.c) * norm(D.*(data.A * u(m+1:m+n))) / (-data.c'*u(m+1:m+n)) / scale;
    else
        unb_res = inf;
    end
        
    if (-data.b'*u(1:m) < 0)
        inf_res = norm(D.*data.b) * norm(E.*(data.A' * u(1:m) + v(m+1:m+n))) / (data.b'*u(1:m)) / scale;
    else
        inf_res = inf;
    end
    
    ratio = mu/eps;
    error_ratio = max(gap,max(err_pri,err_dual))/eps;
    solved = error_ratio < 1;
    infeasible = inf_res < eps;
    unbounded = unb_res < eps;
        
%     ttime = toc;
    
%     fprintf('i: %5d, mu: %3.2e, k: %5d presi: %3.7e dresi: %3.7e, dgap: %3.7e, time: %3.2e \n', ...
%                 i, mu, k, err_pri, err_dual, gap, ttime);
    
    if (solved || infeasible || unbounded)
        break;
    end
    
    if ratio>10
        gamma = 2;
    elseif ratio>1 && ratio<=10
        gamma = 1;
    elseif ratio>0.5 && ratio<=1
        gamma = 0.9;
    elseif ratio>0.1 && ratio<=0.5
        gamma = 0.8;
    elseif ratio>0.05 && ratio<=0.1
        gamma = 0.7;
    elseif ratio>0.01 && ratio<=0.05
        gamma = 0.6;
    elseif ratio>0.005 && ratio<=0.01
        gamma = 0.5;
    elseif ratio>0.001 && ratio<=0.005
        gamma = 0.4;
    else
        gamma = 0.3;
    end
    
    if (k+1>max_iters)
        break;
    end
    
    if error_ratio>6 && error_ratio<=10
        sigma = 0.5;  
    elseif error_ratio>3 && error_ratio<=6
        sigma = 0.6; 
    elseif error_ratio>1 && error_ratio<=3 && ratio > 0.1
        sigma = 0.7;
        gamma = gamma*0.4;
        final_check = 1;
    elseif error_ratio>1 && error_ratio<=3 && ratio < 0.1
        sigma = 0.8;
        gamma = gamma*0.4;
        final_check = 1;
    end
        
    mu = mu*sigma; 
    
    %% reinitialization
    v(m+1:l) = (mu/beta)./u(m+1:l);
    
    %% reinitialize beta
    if adaptive
        beta = 1; 
        v(m+1:l) = (mu/beta)./u(m+1:l);
        u(m+1:l) = (mu/beta)./v(m+1:l);
        beta = BBspectral(work, data, mu, n, m, l, k, u, v, rho_y, use_indirect, extra_verbose, ...
            h, g, gTh, beta, alpha, eps_cor, eps_pen);
        v(m+1:l) = (mu/beta)./u(m+1:l);
    end
    
end

if (k+1 > max_iters); k=k+1; end
if (i+1 == max_outiters); i=i+1; end
ttime = toc;

%% Certificate of infeasibility
tau = abs(u(end));
kap = abs(v(end)) / (sc_b * sc_c * scale);

y = u(1:m) / tau;
x = u(m+1:m+n) / tau;
s = v(m+1:m+n) / tau;

if (tau > undet_tol)
    status = 'solved';
else
    y = nan(m,1);
    x = nan(n,1);
    s = nan(n,1);
    
    y_h = u(1:m);
    x_h = u(m+1:m+n);
    s_h = v(m+1:m+n);
    if norm((u+ut)/2)<=2*(undet_tol*sqrt(l))
        status = 'undetermined';
    elseif data.c'*x_h < data.b'*y_h
        status = 'infeasible';
        y = y_h * scale * sc_b * sc_c /(data.b'*y_h);
        s = s_h * scale * sc_b * sc_c /(data.b'*y_h);
        x = -x_h * scale * sc_b * sc_c /(data.c'*x_h);
    else
        status = 'unbounded';
        y = y_h * scale * sc_b * sc_c /(data.b'*y_h);
        s = s_h * scale * sc_b * sc_c /(data.b'*y_h);
    end
end

info.status    = status;
info.outiter  = i; 
info.iter = k;

info.resPri    = err_pri;
info.resDual   = err_dual;
info.relGap    = gap;
info.time      = ttime; 

if (normalize)
    x = x ./ (E * sc_b);
    y = y ./ (D * sc_c);   
    s = s .* (E / (sc_c * scale));
end

info.pobj    = data.c'*x; 
end
