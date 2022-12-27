function compile_direct(flags, common_abip)

abip_include = strjoin("-I" + [fullfile("include");
    fullfile("linsys");
    fullfile("external");
    fullfile("external", "amd");
    fullfile("external", "ldl")]);

% If use MKL, then the lib_path is your MKL path. For example, in Windows
lib_path = 'C:\Program Files (x86)\Intel\oneAPI\mkl\2022.1.0\lib\intel64';
% or 
% lib_path = 'C:\Program Files (x86)\Intel\oneAPI\mkl\2021.2.0\lib\intel64';

% in MACOS or *nix
% lib_path = "/opt/intel/oneapi/mkl/2021.2.0/lib/";

if exist(lib_path) % #ok
    platform = computer('arch');
    lib_path = "C:\'Program Files (x86)'\Intel\oneAPI\mkl\2022.1.0\lib\intel64";
    lib_path = "-L" + lib_path;
else
    platform = "nomkl";
end % End if

mkl_macro = "-DABIP_PARDISO";
pardiso_src = fullfile("..", "..", "c", "linsys", "abip_pardiso.c");

if platform == "win64"
    fprintf("Linking MKL in Windows \n");
    flags.link = [flags.link, ' -lmkl_intel_ilp64',...
        ' -lmkl_core', ' -lmkl_sequential '];
elseif platform == "maci64"
    fprintf("Linking MKL in MacOS \n");
    flags.link = [flags.link, ' -lmkl_intel_ilp64',...
        ' -lmkl_core', ' -lmkl_sequential '];
elseif platform == "glnxa64"
    fprintf("Linking MKL in Linux \n");
    flags.link = [flags.link, ' -lmkl_intel_ilp64',...
        ' -lmkl_core', ' -lmkl_sequential '];
else
    lib_path = '';
    mkl_macro = '';
    if mexext == "mexw64"
        flags.link = '-lut -lmwblas -lmwlapack';
    else
        flags.link = '-lm -lut -lmwblas -lmwlapack';
    end % End if
    pardiso_src = '';
    fprintf("Unknown platform. Not using MKL \n");
end % End if

if (flags.COMPILE_WITH_OPENMP)
    cmd = sprintf('mex -O %s %s %s %s COMPFLAGS="/openmp \\$COMPFLAGS" CFLAGS="\\$CFLAGS -fopenmp" %s %s', flags.arr, flags.LCFLAG, flags.INCS, abip_include, flags.INT);
else
    cmd = sprintf ('mex -O %s %s %s %s %s %s %s %s', mkl_macro, flags.arr, lib_path, flags.LCFLAG, flags.INCS, abip_include, flags.INT);
end

ldl_path = fullfile("external", "ldl");
amd_path = fullfile("external", "amd");

ldl_files = ["ldl.c"];
amd_files = ["amd_order", "amd_dump", "amd_postorder", "amd_post_tree", ...
    "amd_aat", "amd_2", "amd_1", "amd_defaults", "amd_control", ...
    "amd_info", "amd_valid", "amd_global", "amd_preprocess"];
amd_files = amd_files + ".c";

abip_ldl = fullfile(ldl_path, ldl_files);
abip_ldl = strjoin(abip_ldl);
abip_amd = fullfile(amd_path, amd_files);
abip_amd = strjoin(abip_amd);

cmd = sprintf("%s %s %s %s %s", cmd, abip_amd, abip_ldl,...
    fullfile("linsys", "direct.c"), pardiso_src);
cmd = sprintf ('%s %s %s %s %s -output abip_direct', cmd, common_abip, flags.link, flags.LOCS, flags.BLASLIB);

eval(cmd);

