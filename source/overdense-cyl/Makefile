# This is meant to be included by the osiris Makefile
ifndef DISABLE_OVERDENSE_CYL
# Only proceed if the required sub-modules are enabled
ifndef DISABLE_CYLMODES
ifndef DISABLE_OVERDENSE

FPPF += -DENABLE_OVERDENSE_CYL

# overdense-cyl algorithm source files
src += os-spec-overdense-cyl.f03 os-particles-overdense-cyl.f03 \
	   os-simulation-overdense-cyl.f03


# Add local path to vpath
vpath %.f03 overdense-cyl

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