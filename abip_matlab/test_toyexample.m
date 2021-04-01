data.A = sparse(A);
data.b = full(b);
data.c = full(c);
params_copl = struct('max_outiters', 100, 'max_iters', 10000);
[x, y, s, info_copl] = copl_matlab_large(data,params_copl);
