function test_lasso
clear
addpath('..');

title = ["m","n","seed","lambda","eps","abipsocp time","abipsocp iter","abipsocp obj_p","abipsocp sgm","abipqp time","abipqp iter","abipqp obj_p","abipqp sgm",...
    "scsqp time","scsqp iter","scsqp obj_p","scsqp sgm","scscp time","scscp iter","scscp obj_p","scscp sgm",...
    "grbsocp time","grbsocp iter","grbsocp obj_p","grbsocp sgm","grbqp time","grbqp iter","grbqp obj_p","grbqp sgm"];
format = ["%d","%d","%d","%.1e","%.1e","%f","%d","%f","%f","%f","%d","%f","%f","%f","%d","%f","%f","%f","%d","%f","%f","%f","%d","%f","%f","%f","%d","%f","%f"];
results_file = struct();
count = 1;
eps = 1e-3;
time_limit = 2000;
csv_file = 'LASSO';
csv_path = './';
logfile = ['LASSO_', date(), '_log.txt'];
sh = 10;

abip_time_list = [];
abipqp_time_list = [];
scs_time_list = [];
scscp_time_list = [];
grb_time_list = [];
grbqp_time_list = [];

log_flag = 1;
csv_flag = 1;
seed = 1;

if log_flag
    diary(logfile);
    diary on
end 


for m_lasso = [1000,2000,5000]
    for n_lasso = [5000,10000,15000]
        


        [A,b,lambda] = get_lasso_simu_data(m_lasso,n_lasso,seed);


        abip_error_flag = 0;
        abipqp_error_flag = 0;
        scs_error_flag = 0;
        scscp_error_flag = 0;
        gurobi_error_flag = 0;
        gurobiqp_error_flag = 0;
       
        results_file(count).m = m_lasso;
        results_file(count).n = n_lasso;
        results_file(count).seed = seed;
        results_file(count).lambda  = lambda;
        results_file(count).eps = eps;
        results_file(count).abip_time      = inf;
        results_file(count).abip_iter      = inf;
        results_file(count).abip_obj_p     = inf;
        results_file(count).abip_sgm     = inf;
        results_file(count).abipqp_time      = inf;
        results_file(count).abipqp_iter      = inf;
        results_file(count).abipqp_obj_p     = inf;
        results_file(count).abipqp_sgm     = inf;
        results_file(count).scs_time      = inf;
        results_file(count).scs_iter      = inf;
        results_file(count).scs_obj_p     = inf;
        results_file(count).scs_sgm     = inf;
        results_file(count).scscp_time      = inf;
        results_file(count).scscp_iter      = inf;
        results_file(count).scscp_obj_p     = inf;
        results_file(count).scscp_sgm     = inf;
        results_file(count).gurobi_time    = inf;
        results_file(count).gurobi_iter    = inf;
        results_file(count).gurobi_obj_p   = inf;
        results_file(count).gurobi_sgm   = inf;
        results_file(count).gurobiqp_time    = inf;
        results_file(count).gurobiqp_iter    = inf;
        results_file(count).gurobiqp_obj_p   = inf;
        results_file(count).gurobiqp_sgm   = inf;
    
        
        
%         %ABIP socp
        
        abip_settings = struct('prob_type', 0, 'linsys_solver',1, 'time_limit', time_limit,'eps_p',eps);
        
        data.X = A;
        data.y = b;
        data.lambda = lambda;

        try
            [abip_sol,abip_info] = abip_ml(data,abip_settings);
        catch
            abip_error_flag = 1;
            results_file(count).abip_time      = time_limit;
            results_file(count).abip_iter      = inf;
            results_file(count).abip_obj_p     = inf;
        end
        clear data;
        
        % ABIP QP
        try
            abipqp_data.Q = [speye(m_lasso),sparse(m_lasso,2*n_lasso);sparse(2*n_lasso,m_lasso+2*n_lasso)];
            abipqp_data.c = [zeros(m_lasso,1);lambda*ones(2*n_lasso,1)];
            abipqp_data.b = b;
            abipqp_data.A = [speye(m_lasso),A,-A];
            
            abipqp_cones.f = m_lasso;
            abipqp_cones.l = 2 * n_lasso;
            
            abipqp_settings = struct('linsys_solver',4, 'time_limit', time_limit,'eps_p',eps,'eps_inf',eps*1e-2,'eps_unb',eps*1e-2);
            [abipqp_sol,abipqp_info] = abip_qcp(abipqp_data,abipqp_cones,abipqp_settings);
        catch
            abipqp_error_flag = 1;
            results_file(count).abipqp_time      = time_limit;
            results_file(count).abipqp_iter      = inf;
            results_file(count).abipqp_obj_p     = inf;
        end
        clear abipqp_data;
        clear abipqp_cones;
          

         %scs_socp
         try
            scscp_data.A = [1,sparse(1,3+m_lasso+2*n_lasso);
                            sparse(m_lasso,2),A,-A,sparse(m_lasso,2),speye(m_lasso);
                            1,1,sparse(1,2*n_lasso),-sqrt(2),0,sparse(1,m_lasso);
                            -1,1,sparse(1,2*n_lasso),0,-sqrt(2),sparse(1,m_lasso)];

            scscp_data.A = sparse(scscp_data.A');

            scscp_data.b = [0;1;lambda*ones(2*n_lasso,1);0;0;zeros(m_lasso,1)];

            scscp_data.c = -[1;b;zeros(2,1)];
            scscp_cones.z = 2;
            scscp_cones.l = 2*n_lasso;
            scscp_cones.q = [2+m_lasso];

            scscp_settings = struct('eps_abs', eps, 'eps_rel', eps, 'time_limit_secs', time_limit,'max_iters',10e10);
            [~, ~, ~, res_scscp] = scs(scscp_data, scscp_cones, scscp_settings);
        catch
            scscp_error_flag = 1;
            results_file(count).scscp_time     = time_limit;
            results_file(count).scscp_iter     = inf;
            results_file(count).scscp_pobj    = inf;
        end
        clear scscp_data;
        clear scscp_cones;

        %scs_qp
        try
            scs_data.P = spdiags([ones(m_lasso,1);zeros(2*n_lasso,1)],0,sparse(m_lasso+2*n_lasso,m_lasso+2*n_lasso));
            scs_data.c = [zeros(m_lasso,1);lambda*ones(2*n_lasso,1)];
            scs_data.A = sparse([speye(m_lasso),A,-A;
                            sparse(n_lasso,m_lasso),-speye(n_lasso),sparse(n_lasso,n_lasso);
                            sparse(n_lasso,m_lasso),sparse(n_lasso,n_lasso),-speye(n_lasso)]);
            scs_data.b = [b;zeros(2*n_lasso,1)];
            scs_cones.z = m_lasso;
            scs_cones.l = 2*n_lasso;

            scs_settings = struct('eps_abs', eps, 'eps_rel', eps, 'time_limit_secs', time_limit,'max_iters',10e10);
            [~, ~, ~, res_scs] = scs(scs_data, scs_cones, scs_settings);
        catch
            scs_error_flag = 1;
            results_file(count).scs_time     = time_limit;
            results_file(count).scs_iter     = inf;
            results_file(count).scs_obj_p    = inf;
        end
         clear scs_data;
         clear scs_cones;
         

        %gurobi socp
        try
            A1  = [1, 0        , sparse(1, m_lasso+2*n_lasso);
                  sparse(m_lasso, 2), speye(m_lasso), A, -A  ];
            [m,n] = size(A1);

            nz = 2 * n_lasso;
            nx = n - nz - 2;
            c  = [0; 1 ; sparse(m_lasso, 1); lambda*ones(2*n_lasso, 1)];
            model.obj   = full(c);
            model.A     = sparse(A1);
            model.rhs   = [1; b];
            model.sense = '=';
            model.lb    = [zeros(n-nx-nz,1); -inf*ones(nx,1); zeros(nz,1)];
            quadcon = struct();

            quadcon.Qc = sparse(n, n);
            quadcon.Qc(1,2) = -1;
            quadcon.Qc(2,1) = -1;
            quadcon.Qc(3:2+nx,3:2+nx) = speye(nx);
            quadcon.q = sparse(n, 1);
            quadcon.rhs = 0;

            model.quadcon = quadcon;
            model.modelsense = 'min';

            params.OutputFlag     = 1;
            params.BarConvTol     = eps;
            params.BarQCPConvTol  = eps;
            params.TimeLimit = time_limit;

            res_gurobi = gurobi(model, params);
        catch
            gurobi_error_flag = 1;
            results_file(count).gurobi_time    = time_limit;
            results_file(count).gurobi_iter    = inf;
            results_file(count).gurobi_obj_p   = inf;
        end
        
        clear A1;
        clear model;
        clear quadcon;
        
        %gurobi qp
        
        try
            grbqp_model.Q = 0.5 * spdiags([ones(m_lasso,1);zeros(2*n_lasso,1)],0,sparse(m_lasso+2*n_lasso,m_lasso+2*n_lasso));
            grbqp_model.obj = [zeros(m_lasso,1);lambda*ones(2*n_lasso,1)];
            grbqp_model.A = [speye(m_lasso),A,-A;];
            grbqp_model.sense = repmat('=',m_lasso,1);
            grbqp_model.modelsense = 'min';
            grbqp_model.rhs = b;
            grbqp_model.lb = [-inf*ones(m_lasso,1);zeros(2*n_lasso,1)];
            grbqp_model.ub = inf * ones(m_lasso + 2*n_lasso,1);
            
            grbqp_params.BarConvTol     = eps;
            grbqp_params.BarQCPConvTol  = eps;
            grbqp_params.TimeLimit = time_limit;

            res_gurobiqp = gurobi(grbqp_model, grbqp_params);
        catch
            gurobiqp_error_flag = 1;
            results_file(count).gurobiqp_time    = time_limit;
            results_file(count).gurobiqp_iter    = inf;
            results_file(count).gurobiqp_obj_p   = inf;
        end
        clear grbqp_model;
        
        
        % result
    
        if(abip_error_flag == 0)
            try
                results_file(count).abip_time      = abip_info.runtime;
                results_file(count).abip_iter      = abip_info.admm_iter;
                results_file(count).abip_obj_p     = abip_info.pobj;
            catch
                results_file(count).abip_time      = time_limit;
                results_file(count).abip_iter      = inf;
                results_file(count).abip_obj_p     = inf;
            end
        end
        if(abipqp_error_flag == 0)
            try
                results_file(count).abipqp_time      = abipqp_info.runtime;
                results_file(count).abipqp_iter      = abipqp_info.admm_iter;
                results_file(count).abipqp_obj_p     = abipqp_info.pobj;
            catch
                results_file(count).abipqp_time      = time_limit;
                results_file(count).abipqp_iter      = inf;
                results_file(count).abipqp_obj_p     = inf;
            end
        end
        if(scs_error_flag == 0)
            try
                results_file(count).scs_time     = (res_scs.setup_time + res_scs.solve_time) / 1e3;
                results_file(count).scs_iter     = res_scs.iter;
                results_file(count).scs_obj_p    = res_scs.pobj;
            catch
                results_file(count).scs_time     = time_limit;
                results_file(count).scs_iter     = inf;
                results_file(count).scs_obj_p    = inf;
            end
        end
        if(scscp_error_flag == 0)
            try
                results_file(count).scscp_time     = (res_scscp.setup_time + res_scscp.solve_time) / 1e3;
                results_file(count).scscp_iter     = res_scscp.iter;
                results_file(count).scscp_obj_p    = res_scscp.pobj;
            catch
                results_file(count).scscp_time     = time_limit;
                results_file(count).scscp_iter     = inf;
                results_file(count).scscp_obj_p    = inf;
            end
        end
        if(gurobi_error_flag == 0)
            try
                results_file(count).gurobi_time    = res_gurobi.runtime;
                results_file(count).gurobi_iter    = res_gurobi.baritercount;
                results_file(count).gurobi_obj_p   = res_gurobi.objval;
            catch
                results_file(count).gurobi_time    = time_limit;
                results_file(count).gurobi_iter    = inf;
                results_file(count).gurobi_obj_p   = inf;
            end
        end
        if(gurobiqp_error_flag == 0)
            try
                results_file(count).gurobiqp_time    = res_gurobiqp.runtime;
                results_file(count).gurobiqp_iter    = res_gurobiqp.baritercount;
                results_file(count).gurobiqp_obj_p   = res_gurobiqp.objval;
            catch
                results_file(count).gurobiqp_time    = time_limit;
                results_file(count).gurobiqp_iter    = inf;
                results_file(count).gurobiqp_obj_p   = inf;
            end
        end
        
        
        
        abip_time_list = [abip_time_list, results_file(count).abip_time];
        abipqp_time_list = [abipqp_time_list, results_file(count).abipqp_time];
        scs_time_list = [scs_time_list, results_file(count).scs_time];
        scscp_time_list = [scscp_time_list, results_file(count).scscp_time];
        grb_time_list = [grb_time_list, results_file(count).gurobi_time];
        grbqp_time_list = [grbqp_time_list, results_file(count).gurobiqp_time];

        abip_sgm = calculate_SGM(abip_time_list, sh);
        abipqp_sgm = calculate_SGM(abipqp_time_list, sh);
        scs_sgm = calculate_SGM(scs_time_list, sh);
        scscp_sgm = calculate_SGM(scscp_time_list, sh);
        grb_sgm = calculate_SGM(grb_time_list, sh);
        grbqp_sgm = calculate_SGM(grbqp_time_list, sh);
        
        min_sgm = min([abip_sgm, abipqp_sgm, scs_sgm, scscp_sgm, grb_sgm, grbqp_sgm]);

        results_file(count).abip_sgm  = abip_sgm / min_sgm;
        results_file(count).abipqp_sgm  = abipqp_sgm / min_sgm;
        results_file(count).scs_sgm   = scs_sgm / min_sgm;
        results_file(count).scscp_sgm   = scscp_sgm / min_sgm;
        results_file(count).gurobi_sgm = grb_sgm / min_sgm;
        results_file(count).gurobiqp_sgm = grbqp_sgm / min_sgm;
        
        if csv_flag
            if(count > 1)
                delete(current_csv_file);
            end
            current_csv_file = [csv_path,csv_file,num2str(count),date(),'.csv'];
            write2csv(current_csv_file, title, results_file, format, "struct");
            count = count + 1;
        end
        
    clear A;
    clear b;
    end
end

if log_flag
    diary off
end

if csv_flag
    filename = [csv_path,csv_file,date(),'.csv'];
    write2csv(filename, title, results_file, format, "struct");
end
end
