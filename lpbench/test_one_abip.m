% @date: 06/23/2023
% @author: zcw, gwz
% function to invoke in cmd mode;
%  - fdir : directory of the mps file
%  - fodir: directory to save
%  - eps  : problem accuracy, see the paper for details.
% example run at matlab:
%   test_one_abip("*.mps", '/tmp/', 1e-4, 1000)


function [infolp] = test_one_abip(fdir, fodir, eps, timelimit)
% other method attr comes from the loaded .mat file
% with `params` preset with default settings.
params = abip_get_params();
params.file_dir = fdir;
params.output_dir = fodir;
params.tol = eps;
params.timelimit = timelimit;
K.l = 1;

if isfile(params.file_dir)
  try
    fprintf("=== running  : %s\n", params.file_dir);
    fi = params.file_dir;
    fname_arr = strsplit(params.file_dir, '/');
    fname = fname_arr{1, length(fname_arr)};
    fo = strcat(params.output_dir, '/', fname, sprintf('.%.0e', eps), '.json');

    data = preprocess(fi);
    [xlp, ylp, slp, infolp] = abip(data, K, params);
    infolp.xlp = xlp;
    infolp.ylp = ylp;
    fprintf("===  start saving to: %s\n", fo);
    output_json = jsonencode(infolp);
    file_o = fopen(fo, 'w');
    fprintf(file_o, output_json);
    fprintf("=== finish saving to: %s\n", fo);
    %% add primal dual solutions.
    % save primal
    fprimal_name = strcat(params.output_dir, '/', fname, '.primal.txt');
    fprimal = fopen(fprimal_name, 'w');
    for xi = xlp
      fprintf(fprimal, sprintf('%.4f\n', xi));
    end
    % save dual
    fdual_name = strcat(params.output_dir, '/', fname, '.dual.txt');
    fdual = fopen(fdual_name, 'w');
    for xi = ylp
      fprintf(fdual, sprintf('%.4f\n', xi));
    end
  catch
    warning(sprintf("!!! %s failed", fdir));
  end
else
  % else if an directory to run
  disp('not support a directory to run')
end
end
