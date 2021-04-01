ifeq ($(OS),Windows_NT)
UNAME = CYGWINorMINGWorMSYS
else
UNAME = $(shell uname -s)
endif

CC = gcc
# For cross-compiling with mingw use these.
#CC = i686-w64-mingw32-gcc -m32
#CC = x86_64-w64-mingw32-gcc-4.8

ifneq (, $(findstring CYGWIN, $(UNAME)))
ISWINDOWS := 1
else
ifneq (, $(findstring MINGW, $(UNAME)))
ISWINDOWS := 1
else
ifneq (, $(findstring MSYS, $(UNAME)))
ISWINDOWS := 1
else
ifneq (, $(findstring mingw, $(CC)))
ISWINDOWS := 1
else
ISWINDOWS := 0
endif
endif
endif
endif

ifeq ($(UNAME), Darwin)
# we're on apple, no need to link rt library
LDFLAGS += -lm
SHARED = dylib
SONAME = -install_name
else
ifeq ($(ISWINDOWS), 1)
# we're on windows (cygwin or msys)
LDFLAGS += -lm
SHARED = dll
SONAME = -soname
#TODO: probably doesn't work:
else
# we're on a linux system, use accurate timer provided by clock_gettime()
LDFLAGS += -lm -lrt
SHARED = so
SONAME = -soname
endif
endif

# Add on default CFLAGS
CFLAGS += -g -Wall -Wwrite-strings -pedantic -O3 -funroll-loops -Wstrict-prototypes -I. -Iinclude -Ilinsys
ifneq ($(ISWINDOWS), 1)
CFLAGS += -fPIC
endif

LINSYS = linsys
DIRSRC = $(LINSYS)/direct
DIRSRCEXT = $(DIRSRC)/external

OUT = out
AR = ar
ARFLAGS = rv
ARCHIVE = $(AR) $(ARFLAGS)
RANLIB = ranlib

OPT_FLAGS =
########### OPTIONAL FLAGS ##########
DLONG = 0
ifneq ($(DLONG), 0)
OPT_FLAGS += -DDLONG=$(DLONG) # use longs rather than ints
endif
CTRLC = 1
ifneq ($(CTRLC), 0)
OPT_FLAGS += -DCTRLC=$(CTRLC) # graceful interrupts with ctrl-c
endif
SFLOAT = 0
ifneq ($(SFLOAT), 0)
OPT_FLAGS += -DSFLOAT=$(SFLOAT) # use floats rather than doubles
endif
NOVALIDATE = 0
ifneq ($(NOVALIDATE), 0)
OPT_FLAGS += -DNOVALIDATE=$(NOVALIDATE)$ # remove data validation step
endif
NOTIMER = 0
ifneq ($(NOTIMER), 0)
OPT_FLAGS += -DNOTIMER=$(NOTIMER) # no timing, times reported as nan
endif
COPYAMATRIX = 1
ifneq ($(COPYAMATRIX), 0)
OPT_FLAGS += -DCOPYAMATRIX=$(COPYAMATRIX) # if normalize, copy A
endif

### VERBOSITY LEVELS: 0,1,2
EXTRA_VERBOSE = 0
ifneq ($(EXTRA_VERBOSE), 0)
OPT_FLAGS += -DEXTRA_VERBOSE=$(EXTRA_VERBOSE) # extra verbosity level
endif

MATLAB_MEX_FILE = 0
ifneq ($(MATLAB_MEX_FILE), 0)
OPT_FLAGS += -DMATLAB_MEX_FILE=$(MATLAB_MEX_FILE) # matlab mex
endif
PYTHON = 0
ifneq ($(PYTHON), 0)
OPT_FLAGS += -DPYTHON=$(PYTHON) # python extension
endif

# debug to see var values, e.g. 'make print-OBJECTS' shows OBJECTS value
print-%: ; @echo $*=$($*)
