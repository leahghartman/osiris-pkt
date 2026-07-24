# This is meant to be included by the osiris Makefile
ifdef ENABLE_CUDA

# include new source files for compilation
# we list the cuda files separately for a linking step in main makefile

cuda_src = os-sys-cuda.cu \
	os-vdf-cuda.cu \
	os-chunk-pool-cuda.cu \
	os-mem-transfer-manager-cuda.cu \
	os-species-chunk-manager-cuda.cu \
	os-spec-push-cuda.cu \
	os-nconf-cuda.cu \
	os-current-cuda.cu \
	os-emf-cuda.cu \
	os-species-cuda.cu \
	os-sim-cuda.cu \
	os-cuda-fortran-c-calls.cu \
	os-sort-cuda.cu \
	os-dlb-cuda.cu \
	os-spec-validate-cuda.cu

src += $(cuda_src) \
	os-node-conf-cuda.f03 \
	os-current-define-cuda.f03 \
	os-spec-cuda.f03 \
	os-part-cuda.f03 \
	os-sim-tiles-cuda.f03 \
	os-sim-group-cuda.f03

# Add local path to vpath
vpath %.f03 cuda
vpath %.cu cuda

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif

endif
