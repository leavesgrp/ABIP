function test_cblib
addpath('..');

data_path = 'path to your cblib data path';

files = dir([data_path '/*.gz']);


title  = ["Problem","m","n","eps","abip time","abip iter","abip pobj","abip SGM",...,
          "scs time","scs iter","scs pobj", "scs SGM"];
format = ["%s","%d","%d","%.1e","%.2e","%d","%f","%f","%.2e","%d","%f","%f"];

output_flag = 1;
eps = 1e-4;
time_limit = 100;
iter_limit = 1e10;
rho_x = 1;
count = 1;
results_file = struct();
csv_path = '.';
logfile = './cblib_log.txt';
sh = 10;
abip_time_list = [];
scs_time_list = [];
                                                                                                                                
log_flag = output_flag;
csv_flag = output_flag;

if log_flag
    diary(logfile);
    diary on
end

for j = 1:length(files)
  
    cbf_data = files(j).name;

    mosek_data_str = ['read(', data_path, cbf_data, ')'];
            
    [r,mosek_res] = mosekopt(mosek_data_str);
   
    results_file(count).problem   = cbf_data(1:end-7);
    results_file(count).m         = size(mosek_res.prob.a,1);
    results_file(count).n         = size(mosek_res.prob.a,2);
    results_file(count).eps       = eps;

 
    results_file(count).abip_time     = inf;
    results_file(count).abip_iter     = inf;
    results_file(count).abip_setuptime    = inf;
    results_file(count).abip_time_per_iter     = inf;
    results_file(count).abip_pobj     = inf;
    results_file(count).abip_sgm      = inf;
    abip_error_flag = 0;
    
    results_file(count).scs_time     = inf;
    results_file(count).scs_iter     = inf;
    results_file(count).scs_setuptime    = inf;
    results_file(count).scs_time_per_iter     = inf;
    results_file(count).scs_pobj     = inf;
    results_file(count).scs_sgm      = inf;
    scs_error_flag = 0;
 

    try
        
        [abip_data, abip_cones] = get_abip_data_from_mosek(mosek_res);
        abip_settings = struct('verbose',1,'linsys_solver',1, 'time_limit', time_limit,'eps_p', eps,  'eps_d', eps, 'eps_g', eps, 'eps_inf', 1e-5, 'eps_unb', 1e-5, 'rho_y', 1e-6,'rho_x',rho_x);
        [~,abip_info] = abip_qcp(abip_data, abip_cones, abip_settings);
    
    catch
        abip_error_flag = 1;
        results_file(count).abip_time     = inf;
        results_file(count).abip_iter     = inf;
        results_file(count).abip_setuptime    = inf;
        results_file(count).abip_time_per_iter     = inf;
        results_file(count).abip_pobj     = inf;
        results_file(count).abip_sgm      = inf;
    end
    
    clear abip_data abip_cones;
    
    try
        
        [scs_data, scs_cones] = get_scs_data_from_mosek(mosek_res);
        scs_settings = struct('verbose',0,'check_interval',1,'eps_abs',eps,'eps_rel',eps, 'time_limit_secs',time_limit,'max_iters',iter_limit);
        [~,~,~,scs_info] = scs(scs_data, scs_cones, scs_settings);
    
    catch
        scs_error_flag = 1;
        results_file(count).scs_time     = inf;
        results_file(count).scs_iter     = inf;
        results_file(count).scs_setuptime    = inf;
        results_file(count).scs_time_per_iter     = inf;
        results_file(count).scs_pobj     = inf;
        results_file(count).scs_sgm      = inf;
    end
    
    clear scs_data scs_cones;
    
    clear mosek_res;

    
    if(abip_error_flag == 0)
        try
            results_file(count).abip_time     =  abip_info.runtime;
            results_file(count).abip_iter     =  abip_info.admm_iter;
            results_file(count).abip_setuptime    = abip_info.setup_time;
            results_file(count).abip_time_per_iter     = abip_info.solve_time / abip_info.admm_iter;
            results_file(count).abip_pobj     =  abip_info.pobj;
        catch
            results_file(count).abip_time     = time_limit;
            results_file(count).abip_iter     = inf;
            results_file(count).abip_setuptime    = inf;
            results_file(count).abip_time_per_iter     = inf;
            results_file(count).abip_pobj     = inf;
        end
    end
    
    
    if(scs_error_flag == 0)
        try
            results_file(count).scs_time      = (scs_info.setup_time + scs_info.solve_time) / 1e3;
            results_file(count).scs_iter     =  scs_info.iter;
            results_file(count).scs_setuptime    = scs_info.setup_time / 1e3;
            results_file(count).scs_time_per_iter     = scs_info.solve_time / scs_info.iter / 1e3;
            results_file(count).scs_pobj     =  -scs_info.pobj;
        catch
            results_file(count).scs_time     = time_limit;
            results_file(count).scs_iter     = inf;
            results_file(count).scs_setuptime    = inf;
            results_file(count).scs_time_per_iter     = inf;
            results_file(count).scs_pobj     = inf;
        end
    end
    
  
    
    
    abip_time_list = [abip_time_list, results_file(count).abip_time];
    scs_time_list = [scs_time_list, results_file(count).scs_time];
    
    abip_sgm = calculate_SGM(abip_time_list, sh);
    scs_sgm = calculate_SGM(scs_time_list, sh);
    
    
    results_file(count).abip_sgm  = abip_sgm / min([abip_sgm, scs_sgm]);
    results_file(count).scs_sgm   = scs_sgm / min([abip_sgm, scs_sgm]);
    
    
    if csv_flag
        if(count > 1)
            delete(current_csv_file);
        end

        current_csv_file = [csv_path,'cblib',num2str(rho_x),'_',num2str(count) ,'.csv'];
        write2csv(current_csv_file, title, results_file, format, "struct");
    end
    
    
    count = count + 1;

end

if log_flag
    diary off
end

if csv_flag
    filename = [csv_path,'cblib',num2str(rho_x),'_',date(),'.csv'];
    write2csv(filename, title, results_file, format, "struct");
end
end