# This is meant to be included by the osiris Makefile
ifndef DISABLE_NEUTRAL_SPIN
FPPF += -DENABLE_NEUTRAL_SPIN

# TDSE ionization algorithm source files
src += os-neutral-spin.f03 os-particles-neutral-spin.f03 os-simulation-neutral-spin.f03

# Add local path to vpath
vpath %.f03 neutral_spin

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	@$(MAKE) -C ..

clean:
	@$(MAKE) -C .. clean

endif

endif