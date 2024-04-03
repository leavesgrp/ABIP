%% parameter setting for compiling the ABIP solver. 
gpu = false;    		% compile the gpu version of ABIP. 
float = false;         	% using single precision (rather than double) floating points
int = false;         	% use 32 bit integers for indexing

% WARNING: openmp used in MATLAB can cause errors and crashes, use with caution.  
% openmp parallelizes the matrix multiply for the indirect solver (using CG):
flags.COMPILE_WITH_OPENMP = false;
flags.EXTRA_VERBOSE = 0;

flags.BLASLIB = '';
% MATLAB_MEX_FILE env variable sets blasint to ptrdiff_t. 
flags.LCFLAG  = '-DMATLAB_MEX_FILE -DUSE_LAPACK -DCTRLC=1 -DCOPYAMATRIX'; 
flags.INCS = '';
flags.LOCS = '';


% Common source files
abip_common_src = ["linalg.c"; "adaptive.c"; "cs.c"; "util.c"; "abip.c";
                    "ctrlc.c"; "normalize.c"; "abip_version.c"];
abip_common_linsys = ["common.c"];
abip_mexfile = ["abip_mex.c"];

abip_common_src = fullfile('src', abip_common_src);
abip_common_src = strjoin(abip_common_src);
abip_common_linsys = fullfile('linsys', abip_common_linsys);
abip_mexfile = fullfile('mexfile', abip_mexfile);

common_abip = strcat(abip_common_src, " ", abip_common_linsys, " ", abip_mexfile);

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

delete(fullfile(".", "interface", "*." + mexext));
compile_direct(flags, common_abip);
compile_indirect(flags, common_abip);
movefile(fullfile(".", "*." + mexext), fullfile(".", "interface"));
addpath(fullfile("mexfile"));
addpath(fullfile("interface"));
savepath