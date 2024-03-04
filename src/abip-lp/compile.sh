MKL_HOME=/opt/intel/oneapi/mkl/2021.2.0/lib
OMP_HOME=/opt/intel/oneapi/compiler/2021.2.0/mac/compiler/lib
CC=gcc

# Get the operating system name
OS=$(uname -s)

if [ "$OS" = "Linux" ]; then
    echo "This is a Linux system."
    ${CC} -O2 -o abip-direct *.c -DMAKE_BINARY -Wl,--start-group,${MKL_HOME}/libmkl_core.a ${MKL_HOME}/libmkl_intel_lp64.a ${MKL_HOME}/libmkl_intel_thread.a ${OMP_HOME}/libiomp5.a -Wl,--end-group -lm -ldl -lpthread
    ${CC} -O2 -o abip-indirect *.c -DABIP_INDIRECT_SOLVE -DMAKE_BINARY -Wl,--start-group,${MKL_HOME}/libmkl_core.a ${MKL_HOME}/libmkl_intel_lp64.a ${MKL_HOME}/libmkl_intel_thread.a ${OMP_HOME}/libiomp5.a -Wl,--end-group -lm -ldl -lpthread
elif [ "$OS" = "Darwin" ]; then
    echo "This is a macOS system."
    arch -x86_64 ${CC} -O2 -o abip-direct *.c -DMAKE_BINARY ${MKL_HOME}/libmkl_core.a ${MKL_HOME}/libmkl_intel_lp64.a ${MKL_HOME}/libmkl_intel_thread.a ${OMP_HOME}/libiomp5.a -lm -ldl -lpthread
	arch -x86_64 ${CC} -O2 -o abip-indirect *.c -DABIP_INDIRECT_SOLVE -DMAKE_BINARY ${MKL_HOME}/libmkl_core.a ${MKL_HOME}/libmkl_intel_lp64.a ${MKL_HOME}/libmkl_intel_thread.a ${OMP_HOME}/libiomp5.a -lm -ldl -lpthread
else
    echo "This is an unsupported system."
fi
