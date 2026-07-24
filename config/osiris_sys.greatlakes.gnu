##########################################################################
# System specific configuration for OSIRIS
#   System:    GreatLakes at UMich
#   Compilers: GCC family ( gcc, gfortran )
##########################################################################

##########################################################################
# OSIRIS configuration

# modules than can be turned on and off commenting and uncommenting each line
# DISABLE_XXFEL = 1
# DISABLE_TILES = 1
# DISABLE_SHEAR = 1
# DISABLE_RAD = 1
# DISABLE_RADCYL = 1
# DISABLE_CYLMODES = 1
# DISABLE_QED = 1
# DISABLE_QEDCYL = 1
# DISABLE_PGC = 1
# DISABLE_OVERDENSE = 1
# DISABLE_OVERDENSE_CYL = 1
# DISABLE_NEUTRAL_SPIN = 1
# DISABLE_GR = 1

# Set the following to a system optimized timer, or leave commented to use a default one
TIMER = __POSIX_TIMER__

# SIMD
# Uncomment the following line to use SIMD optimized code.
# You can use SIMD = SSE, AVX or AVX2
#SIMD = SSE
#SIMD = AVX
SIMD = AVX2

# Numeric precision (SINGLE|DOUBLE)
# also update fftw flags at bottom of file
PRECISION = SINGLE
#PRECISION = DOUBLE

##########################################################################
# Compilers

F90 = gfortran
F03 = gfortran -std=gnu -ffree-form

cc  = gcc -fstrict-aliasing
CC  = gcc -fstrict-aliasing

F90FLAGS_all = -cpp -pipe -ffree-line-length-none -fno-range-check

# Fortran preprocessor
FPP = gcc -C -E -x assembler-with-cpp -nostdinc

# This flag supresses some hyper-vigilant compilier warnings that have been deemed harmless. The list of warnings
# can be found is ./source/config.mk.warnings. If you are having some strange issues that you need to debug
# comment out DISABLE_PARANOIA and/or read the ./source/config.warnings file to allow the warnings if you think
# they may help you find the issue.
DISABLE_PARANOIA = YES

##########################################################################
# Fortran flags

# External name mangling
UNDERSCORE = FORTRANSINGLEUNDERSCORE

# Flag to enable compilation of .f03 files (so far only intel compiler requires this)
#F03_EXTENSION_FLAG =

# Enable OpenMP code
F90FLAGS_all += --openmp
FPP += -D_OPENMP

# ------------------------- Compilation Targets -------------------------

# Debug
F90FLAGS_debug      = $(F90FLAGS_all) -g -Og -fbacktrace -fbounds-check \
                      -Wall -fimplicit-none -pedantic \
                      -Wimplicit-interface -Wconversion -Wsurprising \
                      -Wunderflow -ffpe-trap=invalid,zero,overflow \
#                      -fsanitize=address

# Production
F90FLAGS_production = $(F90FLAGS_all) -O3 -march=native -flto

# Profile Shark
# you cannot source profile with -ipo -static (which are turned on by -fast)
# After compilation you must generate the symbol table manually (the system gets confused
# because of the extra preprocessor, the C code does not require this)
# go to the directory where the binary is located and run:
# % dsymutil <binary-name>

F90FLAGS_profile    = -g $(F90FLAGS_production)

##########################################################################################
# C flags

# -wd,981,1418,2259 \
#debug
CFLAGS_debug      = -Og -g -Wall -pedantic -march=native -std=c99 \
#                    -fsanitize=address

# profile
CFLAGS_profile    = -g $(CFLAGS_production)

# production
CFLAGS_production = -O3 -march=native -std=c99 -flto

##########################################################################################
# Linker flags

#LDFLAGS =

##########################################################################################
# Libraries

# MPI
#FPP += -D__HAS_MPI_IN_PLACE__
MPI_FCOMPILEFLAGS = $(shell mpifort -show | cut -f 1 --complement -d ' ')
MPI_FLINKFLAGS    = $(shell mpifort -show | cut -f 1 --complement -d ' ')

# Uncomment the following if HDF5 supports parallel MPIIO
H5_HAVE_PARALLEL = 1
# HDF5
H5_ROOT          = ${HDF5_ROOT}
H5_FCOMPILEFLAGS = -I${HDF5_INCLUDE}
H5_FLINKFLAGS    = -L${HDF5_LIB} \
                   -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 \
                   -lz -lrt -ldl -lm -Wl,-rpath -Wl,${HDF5_LIB}

FFTW_FCOMPILEFLAGS = -I${FFTW_INCLUDE}
# flags for single precision
FFTW_FLINKFLAGS    = -L${FFTW_LIB} \
                     -lfftw3f -lm
# flags for double precision
#FFTW_FLINKFLAGS   = -L${FFTW_LIB} \
#                    -lfftw3 -lm

INCPATH = ${MPI_FCOMPILEFLAGS}
