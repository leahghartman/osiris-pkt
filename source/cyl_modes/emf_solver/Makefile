# This is meant to be included by the osiris Makefile

# emf object source files
src +=  os-emf-solver-cyl-modes-yee.f03 os-emf-solver-cyl-modes-fei.f03 \
		os-emf-solver-cyl-modes.f03

# Add local path to vpath
vpath %.f03 cyl_modes/emf_solver

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
