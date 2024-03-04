% @date: 06/23/2023
% @author: zcw, gwz
% function to invoke in cmd mode;
%  - fdir : directory of the mps file
%  - fodir: directory to save
%  - eps  : problem accuracy, see the paper for details.
% example run at matlab:
%   test_one_abip("*.mps", '/tmp/', 1e-4, 1000)


function [data] = save_abip_mps(fdir, fodir)
% other method attr comes from the loaded .mat file
% with `params` preset with default settings.
% params = abip_get_params();
params.file_dir = fdir;
params.output_dir = fodir;
K.l = 1;

if isfile(params.file_dir)
  try
    fprintf("=== running  : %s\n", params.file_dir);
    fi = params.file_dir;
    fname_arr = strsplit(params.file_dir, '/');
    fname = fname_arr{1, length(fname_arr)};
    fo = strcat(params.output_dir, '/', fname);

    data = preprocess(fi);
    data.obj = data.c;
    data.rhs = data.b;
    data.sens = '=';
    if isnan(data.objcon)
        data.objcon = 0.0
    end
    gurobi_write(data, fo)
  catch
    warning(sprintf("!!! %s failed", fdir));
  end
else
  % else if an directory to run
  disp('not support a directory to run')
end
end
