# This is meant to be included by the osiris Makefile
ifndef DISABLE_RADCYL
# Only proceed if the required sub-modules are enabled
ifndef DISABLE_CYLMODES
ifndef DISABLE_RAD

FPPF += -DENABLE_RADCYL

# qed-cyl algorithm source files
src += os-species-radio-cyl.f03 os-particles-radio-cyl.f03 os-simulation-radio-cyl.f03 

# Add local path to vpath
vpath %.f03 radio-cyl

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	@$(MAKE) -C ..

clean:
	@$(MAKE) -C .. clean

endif

endif
endif
endif
