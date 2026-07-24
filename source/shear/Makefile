# This is meant to be included by the osiris Makefile
ifndef DISABLE_SHEAR
FPPF += -DENABLE_SHEAR

# SHEAR algorithm source files
src += os-emf-shear.f03	 os-spec-shear.f03 \
			 os-particles-shear.f03 os-simulation-shear.f03

# Add local path to vpath
vpath %.f03 shear

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif

endif