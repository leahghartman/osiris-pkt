# This is meant to be included by the osiris/spec Makefile

# partticle source object source files
src +=  os-psource-std-cyl-modes.f03 os-psource-beam-cyl-modes.f03 \
		os-psource-py-cyl-modes.f03 os-psource-cyl-modes.f03

# Add local path to vpath
vpath %.f03 cyl_modes/psource

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
