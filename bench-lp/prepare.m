clear; clc; close all;


fc = 'datasets/mittelman'
fcpre = 'datasets/s_mittelman'

files = dir(fullfile(".", fc, '/*.mps'));
nfiles = length(files);
fnames = {files.name}';

for i = 1:nfiles
    fname = fnames{i};
    
    try
        data = preprocess(fullfile('.', fc, fname));
        model.A = data.A;
        model.rhs = data.b;
        model.obj = data.c;
        model.sense = '=';
        fprename = fullfile('.', fcpre, ['s_' fname(1:end-3) 'mps'])
        gurobi_write(model, fprename);
        fprintf("%s succeeded. \n", fprename);
    catch
        fprintf("%s failed. \n", fname);
    end % End try
    
end % End for