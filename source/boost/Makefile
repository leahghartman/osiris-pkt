# This is meant to be included by the osiris Makefile

# emf object source files
src += os-boosted-grid-mgr.f03  os-boosted-diag.f03


# Add local path to vpath
vpath %.f03 boost
vpath %.f90 boost

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
