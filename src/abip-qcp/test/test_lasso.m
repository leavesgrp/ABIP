function test_lasso
clear
addpath('..');

title = ["m","n","lambda","eps","abip time","abip iter","abip obj_p","abip sgm",...
    "scs time","scs iter","scs obj_p","scs sgm","gurobi time","gurobi iter","gurobi obj_p","gurobi sgm"];
format = ["%d","%d","%.1e","%.1e","%f","%d","%f","%f","%f","%d","%f","%f","%f","%d","%f","%f"];
results_file = struct();
count = 1;
eps = 1e-3;
time_limit = 2e3;
csv_path = './';
logfile = ['.LASSO_', date(), '_log.txt'];
sh = 10;

abip_time_list = [];
scs_time_list = [];
grb_time_list = [];

log_flag = 1;
csv_flag = 1;

if log_flag
    diary(logfile);
    diary on
end

for m_lasso = [2000,2500,3000,4000,5000,10000]
    for n_lasso = [50,100,150,200,10000,12500,15000,20000]
    

        abip_error_flag = 0;
        scs_error_flag = 0;
        gurobi_error_flag = 0;
        
        A = sprandn(m_lasso ,n_lasso, 0.15);
        v = zeros(n_lasso,1);
        for i = 1:n_lasso
            if rand() > 0.5
                v(i) = normrnd(0,1/n_lasso);
            end
        end
        epsilon = randn(m_lasso,1);
        b = A*v + epsilon;
        lambda_max = norm(A'*b,Inf);
        lambda = lambda_max / 5;
        
        results_file(count).m = m_lasso;
        results_file(count).n = n_lasso;
        results_file(count).lambda  = lambda;
        results_file(count).eps = eps;
        results_file(count).abip_time      = inf;
        results_file(count).abip_iter      = inf;
        results_file(count).abip_obj_p     = inf;
        results_file(count).abip_sgm     = inf;
        results_file(count).scs_time      = inf;
        results_file(count).scs_iter      = inf;
        results_file(count).scs_obj_p     = inf;
        results_file(count).scs_sgm     = inf;
        results_file(count).gurobi_time    = inf;
        results_file(count).gurobi_iter    = inf;
        results_file(count).gurobi_obj_p   = inf;
        results_file(count).gurobi_sgm   = inf;

        
        
        
        
        
        %ABIP
        
        abip_settings = struct('prob_type', 0, 'time_limit', time_limit,'eps',eps);
        
        data.X = A;
        data.y = b;
        data.lambda = lambda;

        try
            [abip_sol,abip_info] = abip_ml(data,abip_settings);
        catch
            abip_error_flag = 1;
            results_file(count).abip_time      = inf;
            results_file(count).abip_iter      = inf;
            results_file(count).abip_obj_p     = inf;
        end

        
        %scs
        try
            scs_data.A = [1,sparse(1,3+m_lasso+2*n_lasso);
                                sparse(m_lasso,2),A,-A,sparse(m_lasso,2),speye(m_lasso);
                                1,1,sparse(1,2*n_lasso),-sqrt(2),0,sparse(1,m_lasso);
                                -1,1,sparse(1,2*n_lasso),0,-sqrt(2),sparse(1,m_lasso)];
                
            scs_data.A = sparse(scs_data.A');

            scs_data.b = [0;1;lambda*ones(2*n_lasso,1);0;0;zeros(m_lasso,1)];

            scs_data.c = -[1;b;zeros(2,1)];
            scs_cones.z = 2;
            scs_cones.l = 2*n_lasso;
            scs_cones.q = [2+m_lasso];

            scs_settings = struct('eps_abs', eps, 'eps_rel', eps, 'time_limit_secs', 2e3);
            [~, ~, ~, res_scs] = scs(scs_data, scs_cones, scs_settings);
        catch
            scs_error_flag = 1;
            results_file(count).scs_time     = inf;
            results_file(count).scs_iter     = inf;
            results_file(count).scs_obj_p    = inf;
        end
         clear scs_data;
        
        %gurobi
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
            params.TimeLimit = 2000;

            res_gurobi = gurobi(model, params);
        catch
            gurobi_error_flag = 1;
            results_file(count).gurobi_time    = inf;
            results_file(count).gurobi_iter    = inf;
            results_file(count).gurobi_obj_p   = inf;
        end
        clear A1;
        clear model;
        clear quadcon;
        
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
        
        
        abip_time_list = [abip_time_list, results_file(count).abip_time];
        scs_time_list = [scs_time_list, results_file(count).scs_time];
        grb_time_list = [grb_time_list, results_file(count).gurobi_time];

        abip_sgm = calculate_SGM(abip_time_list, sh);
        scs_sgm = calculate_SGM(scs_time_list, sh);
        grb_sgm = calculate_SGM(grb_time_list, sh);
        
        min_sgm = min([abip_sgm, scs_sgm, grb_sgm]);

        results_file(count).abip_sgm  = abip_sgm / min_sgm;
        results_file(count).scs_sgm   = scs_sgm / min_sgm;
        results_file(count).gurobi_sgm = grb_sgm / min_sgm;
        
        

        if csv_flag
            if(count > 1)
                delete(current_csv_file);
            end
            current_csv_file = [csv_path,'LASSO_',num2str(count),'.csv'];
            write2csv(current_csv_file, title, results_file, format, "struct");
            count = count + 1;
        end
        
    clear X;
    clear y;
    end
end

if log_flag
    diary off
end

if csv_flag
    filename = [csv_path,'LASSO_',date(),'.csv'];
    write2csv(filename, title, results_file, format, "struct");
end
end