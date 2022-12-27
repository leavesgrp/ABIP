function test_svm
clear

addpath('..');

data_path = "path to your SVM dataset";
addpath(data_path);
data_dir = data_path;

data_set = dir([data_dir '/*']);

title = ["data","m","n","sparsity","lambda","eps","abip time","abip iter","abip obj_p","abip sgm",...
    "scs time","scs iter","scs obj_p","scs sgm","gurobi time","gurobi iter","gurobi obj_p","gurobi sgm"];
format = ["%s","%d","%d","%f","%.1e","%.1e","%f","%d","%f","%f","%f","%d","%f","%f","%f","%d","%f","%f"];
results_file = struct();
count = 1;
eps = 1e-3;
time_limit = 3e3;
csv_path = '.';
logfile = ['./SVM', date(), '_log.txt'];
sh = 10;
lambda = 1e-3;

abip_time_list = [];
scs_time_list = [];
grb_time_list = [];

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
    
    abip_error_flag = 0;
    scs_error_flag = 0;
    gurobi_error_flag = 0;

    results_file(count).data_name = data_name;
    results_file(count).m = N;
    results_file(count).n = n;
    results_file(count).sparsity = sparsity;
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


    C = 1 / (N*lambda);

    %ABIP

    abip_settings = struct('verbose',1,'prob_type', 1, 'linsys_solver',linsys_solver, 'time_limit', time_limit, 'outer_check_period', ocp);

    data.X = X;
    data.y = y;
    data.lambda = C;

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
        scs_settings = struct('eps_abs', eps, 'eps_rel', eps, 'time_limit_secs', time_limit, 'max_iters',10e10);

        scs_data.A = sparse([-1,-1,0,sparse(1,N),sparse(1,N),sqrt(2),0,sparse(1,n);
                              1,-1,0,sparse(1,N),sparse(1,N),0,sqrt(2),sparse(1,n);
                              0,1,0,sparse(1,N),sparse(1,N),0,0,sparse(1,n);
                              sparse(N,1),sparse(N,1),y,speye(N),-speye(N),sparse(N,1),sparse(N,1),spdiags(y,0,sparse(N,N))*X]);

        scs_data.c = -[0;0;1;ones(N,1)];

        scs_data.b = 1/lambda*[lambda;0;0;1/N*ones(N,1);zeros(N,1);0;0;zeros(n,1)];
        scs_cones.z = 3;
        scs_cones.l = 2*N;
        scs_cones.q = 2+n;
        [scs_sol.x, scs_sol.y, scs_sol.s, res_scs] = scs(scs_data, scs_cones, scs_settings);
    catch
        scs_error_flag = 1;
        results_file(count).scs_time      = inf;
        results_file(count).scs_iter      = inf;
        results_file(count).scs_obj_p     = inf;
    end
    clear scs_data
    clear scs_cones
    clear abip_socp_data


    %gurobi
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
        res_gurobi = gurobi(model, params);
    catch
        gurobi_error_flag = 1;
        results_file(count).gurobi_time    = inf;
        results_file(count).gurobi_iter    = inf;
        results_file(count).gurobi_obj_p   = inf;
    end

    clear model;


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