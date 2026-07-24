# This is meant to be included by the osiris Makefile
ifndef DISABLE_OVERDENSE
FPPF += -DENABLE_OVERDENSE

# Damping algorithm source files
src += os-spec-fluid-moments.f03 os-spec-damp.f03 os-spec-splitting.f03 \
       os-spec-overdense.f03 os-particles-overdense.f03 os-simulation-overdense.f03

# Add local path to vpath
vpath %.f03 overdense

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	@$(MAKE) -C ..

clean:
	@$(MAKE) -C .. clean

endif

endif