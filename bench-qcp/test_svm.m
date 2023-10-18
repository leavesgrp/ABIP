function test_svm
clear

addpath('..');

data_path = "path to your SVM dataset";


addpath(data_path);
data_dir = data_path;

data_set = dir([data_dir '/*']);

title = ["data","m","n","sparsity","lambda","eps","abipsocp time","abipsocp iter","abipsocp obj_p","abipsocp sgm", "abipqp time","abipqp iter","abipqp obj_p","abipqp sgm",...
    "scssocp time","scssocp iter","scssocp obj_p","scssocp sgm","scsqp time","scsqp iter","scsqp obj_p","scsqp sgm","gurobisocp time","gurobisocp iter","gurobisocp obj_p","gurobisocp sgm",...
    "gurobiqp time","gurobiqp iter","gurobiqp obj_p","gurobiqp sgm"];
format = ["%s","%d","%d","%f","%.1e","%.1e","%f","%d","%f","%f","%f","%d","%f","%f","%f","%d","%f","%f", "%f","%d","%f","%f","%f","%d","%f","%f","%f","%d","%f","%f"];
results_file = struct();
count = 1;
eps = 1e-3;
time_limit = 3e3;
linsys_solver = 4;
csv_path = './';
logfile = ['SVM', date(), '_log.txt'];
sh = 10;
lambda = 1e-3;
rho_x = 0.1;

abipsocp_time_list = [];
abipqp_time_list = [];
scssocp_time_list = [];
scsqp_time_list = [];
grbsocp_time_list = [];
grbqp_time_list = [];

log_flag = 1;
csv_flag = 1;

if log_flag
    diary(logfile);
    diary on
end

for j = 1:length(data_set)
    
    data_name = data_set(j).name;
    fprintf("loading %s\n",data_name);
    [y,X] = libsvmread(fullfile(data_dir,data_name));
    y(find(y==2)) = -1;

    ind1 = find(sum(X.^2,1)==0);
    X(:,ind1) = [];
    [N,n] = size(X);
    sparsity = nnz(X)/N/n;
    
    abipsocp_error_flag = 0;
    abipqp_error_flag = 0;
    scssocp_error_flag = 0;
    scsqp_error_flag = 0;
    gurobisocp_error_flag = 0;
    gurobiqp_error_flag = 0;

    results_file(count).data_name = data_name;
    results_file(count).m = N;
    results_file(count).n = n;
    results_file(count).sparsity = sparsity;
    results_file(count).lambda  = lambda;
    results_file(count).eps = eps;
    results_file(count).abipsocp_time      = inf;
    results_file(count).abipsocp_iter      = inf;
    results_file(count).abipsocp_obj_p     = inf;
    results_file(count).abipsocp_sgm     = inf;
    results_file(count).abipqp_time      = inf;
    results_file(count).abipqp_iter      = inf;
    results_file(count).abipqp_obj_p     = inf;
    results_file(count).abipqp_sgm     = inf;
    results_file(count).scssocp_time      = inf;
    results_file(count).scssocp_iter      = inf;
    results_file(count).scssocp_obj_p     = inf;
    results_file(count).scssocp_sgm     = inf;
    results_file(count).scsqp_time      = inf;
    results_file(count).scsqp_iter      = inf;
    results_file(count).scsqp_obj_p     = inf;
    results_file(count).scsqp_sgm     = inf;
    results_file(count).gurobisocp_time    = inf;
    results_file(count).gurobisocp_iter    = inf;
    results_file(count).gurobisocp_obj_p   = inf;
    results_file(count).gurobisocp_sgm   = inf;
    results_file(count).gurobiqp_time    = inf;
    results_file(count).gurobiqp_iter    = inf;
    results_file(count).gurobiqp_obj_p   = inf;
    results_file(count).gurobiqp_sgm   = inf;


    C = 1 / (N*lambda);

    %abip socp
    abipsocp_settings = struct('verbose',1,'prob_type', 1, 'linsys_solver',linsys_solver, 'time_limit', time_limit, 'outer_check_period', 2);

    data.X = X;
    data.y = y;
    data.lambda = C;

    try
        [abipsocp_sol,abipsocp_info] = abip_ml(data,abipsocp_settings);
    catch
        abipsocp_error_flag = 1;
        results_file(count).abipsocp_time      = inf;
        results_file(count).abipsocp_iter      = inf;
        results_file(count).abipsocp_obj_p     = inf;
    end
    clear data
    
    %abip qp
    abipqp_settings = struct('rho_y',1e-6, 'rho_x',rho_x,'eps',eps,'prob_type', 3, 'linsys_solver',linsys_solver, 'time_limit', time_limit, 'outer_check_period', 1);

    data.X = X;
    data.y = y;
    data.lambda = lambda;
    try
         [abipqp_results, abipqp_info] = abip_ml(data,abipqp_settings);
    catch
        abipqp_error_flag = 1;
        results_file(count).abipqp_time      = inf;
        results_file(count).abipqp_iter      = inf;
        results_file(count).abipqp_obj_p     = inf;
    end
    clear data

    %scs socp
    try
        scssocp_settings = struct('eps_abs', eps, 'eps_rel', eps, 'time_limit_secs', time_limit, 'max_iters',10e10);

        scssocp_data.A = sparse([-1,-1,0,sparse(1,N),sparse(1,N),sqrt(2),0,sparse(1,n);
                              1,-1,0,sparse(1,N),sparse(1,N),0,sqrt(2),sparse(1,n);
                              0,1,0,sparse(1,N),sparse(1,N),0,0,sparse(1,n);
                              sparse(N,1),sparse(N,1),y,speye(N),-speye(N),sparse(N,1),sparse(N,1),spdiags(y,0,sparse(N,N))*X]);
        scssocp_data.A = scssocp_data.A';
        scssocp_data.c = -[0;0;1;ones(N,1)];

        scssocp_data.b = 1/lambda*[lambda;0;0;1/N*ones(N,1);zeros(N,1);0;0;zeros(n,1)];
        scssocp_cones.z = 3;
        scssocp_cones.l = 2*N;
        scssocp_cones.q = 2+n;
        [~, ~, ~, res_scssocp] = scs(scssocp_data, scssocp_cones, scssocp_settings);
    catch
        scssocp_error_flag = 1;
        results_file(count).scssocp_time      = inf;
        results_file(count).scssocp_iter      = inf;
        results_file(count).scssocp_obj_p     = inf;
    end
    clear scssocp_data
    clear scssocp_cones


    %scsqp qp
    try
        l = n + N + 1;
        scsqp_data.P = 1/lambda*spdiags([lambda*ones(n,1);zeros(N+1,1)],0,sparse(l,l));
        scsqp_data.c = 1/lambda*[zeros(n,1);1/N*ones(N,1);0];
        scsqp_data.A = -sparse([spdiags(y,0,sparse(N,N))*X,speye(N,N),y;sparse(N,n),speye(N),sparse(N,1)]);
        scsqp_data.b = -[ones(N,1);zeros(N,1)];
        scsqp_cones.l = 2*N;
        scsqp_settings = struct('use_indirect',0,'eps_abs', eps, 'eps_rel', eps, 'time_limit_secs', time_limit, 'max_iters',10e10);
        [~, ~, ~, res_scsqp] = scs(scsqp_data, scsqp_cones, scsqp_settings);
    catch
        scsqp_error_flag = 1;
        results_file(count).scsqp_time      = inf;
        results_file(count).scsqp_iter      = inf;
        results_file(count).scsqp_obj_p     = inf;
    end
    clear scsqp_data
    clear scsqp_cones

    %gurobi socp
    try
        model.obj   = [zeros(n,1);0;C*ones(N,1);zeros(N,1);1];
        model.A     = sparse([spdiags(y,0,sparse(N,N))*X,y,speye(N),-speye(N),sparse(N,1)]);
        model.rhs   = ones(N,1);
        model.sense = '=';
        model.lb = [-inf*ones(n+1,1);zeros(2*N,1);-inf];
        model.modelsense = 'min';
        model.quadcon.Qc = 1/2*spdiags([ones(n,1);zeros(2*N+2,1)],0,sparse(2*N+n+2,2*N+n+2));
        model.quadcon.q = sparse([sparse(2*N+n+1,1);-1]);
        model.quadcon.rhs = 0;

        params.OutputFlag     = 1;
        params.QCPDual        = 1;
        params.BarConvTol     = eps;
        params.BarQCPConvTol  = eps;
        params.TimeLimit = time_limit;
        res_gurobisocp = gurobi(model, params);
    catch
        gurobisocp_error_flag = 1;
        results_file(count).gurobisocp_time    = inf;
        results_file(count).gurobisocp_iter    = inf;
        results_file(count).gurobisocp_obj_p   = inf;
    end

    clear model;
    clear params;

    %gurobi qp
    try
        model.Q = 1/2*spdiags([ones(n,1);zeros(2*N+1,1)],0,sparse(2*N+n+1,2*N+n+1));
        model.obj   = [zeros(n,1);0;C*ones(N,1);zeros(N,1)];
        model.A     = sparse([spdiags(y,0,sparse(N,N))*X,y,speye(N),-speye(N)]);
        model.rhs   = ones(N,1);
        model.sense = '=';
        model.lb = [-inf*ones(n+1,1);zeros(2*N,1)];
        model.modelsense = 'min';

        params.OutputFlag     = 1;
        params.QCPDual        = 1;
        params.BarConvTol     = eps;
        params.BarQCPConvTol  = eps;
        params.TimeLimit = time_limit;
        res_gurobiqp = gurobi(model, params);
    catch
        gurobiqp_error_flag = 1;
        results_file(count).gurobiqp_time    = inf;
        results_file(count).gurobiqp_iter    = inf;
        results_file(count).gurobiqp_obj_p   = inf;
    end

    clear model;
    clear params;

    % result
    if(abipsocp_error_flag == 0)
        try
            results_file(count).abipsocp_time      = abipsocp_info.runtime;
            results_file(count).abipsocp_iter      = abipsocp_info.admm_iter;
            results_file(count).abipsocp_obj_p     = abipsocp_info.pobj;
        catch
            results_file(count).abipsocp_time      = time_limit;
            results_file(count).abipsocp_iter      = inf;
            results_file(count).abipsocp_obj_p     = inf;
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
    if(scssocp_error_flag == 0)
        try
            results_file(count).scssocp_time     = (res_scssocp.setup_time + res_scssocp.solve_time) / 1e3;
            results_file(count).scssocp_iter     = res_scssocp.iter;
            results_file(count).scssocp_obj_p    = res_scssocp.pobj;
        catch
            results_file(count).scssocp_time     = time_limit;
            results_file(count).scssocp_iter     = inf;
            results_file(count).scssocp_obj_p    = inf;
        end
    end
    if(scsqp_error_flag == 0)
        try
            results_file(count).scsqp_time     = (res_scsqp.setup_time + res_scsqp.solve_time) / 1e3;
            results_file(count).scsqp_iter     = res_scsqp.iter;
            results_file(count).scsqp_obj_p    = res_scsqp.pobj;
        catch
            results_file(count).scsqp_time     = time_limit;
            results_file(count).scsqp_iter     = inf;
            results_file(count).scsqp_obj_p    = inf;
        end
    end
    if(gurobisocp_error_flag == 0)
        try
            results_file(count).gurobisocp_time    = res_gurobisocp.runtime;
            results_file(count).gurobisocp_iter    = res_gurobisocp.baritercount;
            results_file(count).gurobisocp_obj_p   = res_gurobisocp.objval;
        catch
            results_file(count).gurobisocp_time    = time_limit;
            results_file(count).gurobisocp_iter    = inf;
            results_file(count).gurobisocp_obj_p   = inf;
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


    abipsocp_time_list = [abipsocp_time_list, results_file(count).abipsocp_time];
    abipqp_time_list = [abipqp_time_list, results_file(count).abipqp_time];
    scssocp_time_list = [scssocp_time_list, results_file(count).scssocp_time];
    scsqp_time_list = [scsqp_time_list, results_file(count).scsqp_time];
    grbsocp_time_list = [grbsocp_time_list, results_file(count).gurobisocp_time];
    grbqp_time_list = [grbqp_time_list, results_file(count).gurobiqp_time];

    abipsocp_sgm = calculate_SGM(abipsocp_time_list, sh);
    abipqp_sgm = calculate_SGM(abipqp_time_list, sh);
    scssocp_sgm = calculate_SGM(scssocp_time_list, sh);
    scsqp_sgm = calculate_SGM(scsqp_time_list, sh);
    grbsocp_sgm = calculate_SGM(grbsocp_time_list, sh);
    grbqp_sgm = calculate_SGM(grbqp_time_list, sh);

    min_sgm = min([abipsocp_sgm, scssocp_sgm, grbsocp_sgm, abipqp_sgm, scsqp_sgm, grbqp_sgm]);

    results_file(count).abipsocp_sgm  = abipsocp_sgm / min_sgm;
    results_file(count).scssocp_sgm   = scssocp_sgm / min_sgm;
    results_file(count).gurobisocp_sgm = grbsocp_sgm / min_sgm;
    results_file(count).abipqp_sgm  = abipqp_sgm / min_sgm;
    results_file(count).scsqp_sgm   = scsqp_sgm / min_sgm;
    results_file(count).gurobiqp_sgm = grbqp_sgm / min_sgm;

    if csv_flag
        if(count > 1)
            delete(current_csv_file);
        end
        current_csv_file = [csv_path,'SVM_',num2str(count),'.csv'];
        write2csv(current_csv_file, title, results_file, format, "struct");
        count = count + 1;
    end
        
    clear X;
    clear y;
end

if log_flag
    diary off
end

if csv_flag
    filename = [csv_path,'SVM_',date(),'.csv'];
    write2csv(filename, title, results_file, format, "struct");
end
end