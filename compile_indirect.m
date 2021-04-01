function compile_indirect(flags, common_abip)

% compile indirect
if (flags.COMPILE_WITH_OPENMP) 
    cmd = sprintf('mex -O %s %s %s %s COMPFLAGS="/openmp \\$COMPFLAGS" CFLAGS="\\$CFLAGS -fopenmp" abip/linsys/indirect/private.c %s -Iabip -Iabip/include -Iabip/linsys %s %s %s -output abip_indirect', flags.arr, flags.LCFLAG, common_abip, flags.INCS, flags.link, flags.LOCS, flags.BLASLIB, flags.INT);
else
    cmd = sprintf('mex -O %s %s %s %s abip/linsys/indirect/private.c %s -Iabip -Iabip/include -Iabip/linsys %s %s %s -output abip_indirect',  flags.arr, flags.LCFLAG, common_abip, flags.INCS, flags.link, flags.LOCS, flags.BLASLIB, flags.INT);
end
eval(cmd);
