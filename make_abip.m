%% parameter setting for compiling the ABIP solver. 
gpu = false;    		% compile the gpu version of ABIP. 
float = false;         	% using single precision (rather than double) floating points
int = false;         	% use 32 bit integers for indexing

% WARNING: openmp used in MATLAB can cause errors and crashes, use with caution.  
% openmp parallelizes the matrix multiply for the indirect solver (using CG):
flags.COMPILE_WITH_OPENMP = false;
flags.EXTRA_VERBOSE = 0;

flags.BLASLIB = '-lmwblas -lmwlapack';
flags.LCFLAG  = '-DMATLAB_MEX_FILE -DUSE_LAPACK -DCTRLC=1 -DCOPYAMATRIX'; % MATLAB_MEX_FILE env variable sets blasint to ptrdiff_t. 
flags.INCS = '';
flags.LOCS = '';

common_abip = 'abip/src/linalg.c abip/src/adaptive.c abip/src/cs.c abip/src/util.c abip/src/abip.c abip/src/ctrlc.c abip/src/normalize.c abip/src/abip_version.c abip/linsys/common.c abip_mex.c';

if (~isempty (strfind (computer, '64')))
    flags.arr = '-largeArrayDims';
else
    flags.arr = '';
end

if (isunix && ~ismac)
    flags.link = '-lm -lut -lrt';
elseif  (ismac)
    flags.link = '-lm -lut';
else
    flags.link = '-lut';
    flags.LCFLAG = sprintf('-DNOBLASSUFFIX %s', flags.LCFLAG);
end

if (float)
    flags.LCFLAG = sprintf('-DSFLOAT %s', flags.LCFLAG);
end

if (int)
    flags.INT = '';
else
    flags.INT = '-DDLONG';
end

if (flags.COMPILE_WITH_OPENMP)
    flags.link = strcat(flags.link, ' -lgomp');
end

if (flags.EXTRA_VERBOSE)
    flags.LCFLAG = sprintf('-DEXTRA_VERBOSE %s', flags.LCFLAG);
end

compile_direct(flags, common_abip);
compile_indirect(flags, common_abip);
% if (gpu)
%    compile_gpu(flags, common_abip);
% end

% compile abip_version
mex -O -Iabip/include abip/src/abip_version.c abip_version_mex.c -output abip_version

%%
clear data
load('../netlib/feasible/ADLITTLE.mat');
disp('Example run:');
A = Problem.A; 
b = Problem.b; 
c = Problem.aux.c; 
lbounds = Problem.aux.lo; 
ubounds = Problem.aux.hi; 
[m, n] = size(A);

[A,b,c,info] = presolve(A, b, c, lbounds, ubounds);
data.A = sparse(A);
data.b = full(b);
data.c = full(c);
params_abip = struct('max_ipm_iters', 100, 'max_admm_iters', 1000000, 'adaptive', 1, 'verbose', 1);
[x, y, s, info_abip] = abip_direct(data, params_abip);
[x_indirect, y_indirect, s_indirect, info_abip_indirect] = abip_indirect(data, params_abip);
% 
% if (gpu)
%    [x,y,s,info] = abip_gpu(data,cones,[]);
% end
% 
% test-warm start with solution
% disp('Warm-starting:')
% data.x = x;
% data.y = y;
% data.s = s;
% [x, y, s, info] = abip_indirect(data, params_abip);

disp('SUCCESSFULLY INSTALLED ABIP')
