# This is meant to be included by the osiris Makefile
ifndef DISABLE_XXFEL
FPPF += -DENABLE_XXFEL

# Damping algorithm source files
src += os-simulation-xxfel.f03

# Add local path to vpath
vpath %.f03 xxfel

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	@$(MAKE) -C ..

clean:
	@$(MAKE) -C .. clean

endif

endif