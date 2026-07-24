# This is meant to be included by the osiris Makefile
ifndef DISABLE_CYLMODES
FPPF += -DENABLE_CYLMODES

# Cylindircal Modes algorithm source files
src +=  os-cyl-modes.f03 os-grid-cyl-modes.f03 os-current-cyl-modes.f03 \
		os-emf-bound-cyl-modes.f03

# emf solver object source files
include cyl_modes/emf_solver/Makefile

# boosted frame source files
src += os-boosted-diag-cyl-modes.f03

src +=  os-emf-cyl-modes.f03 os-emf-interp-cm.f03 \
		os-spec-cyl-modes-define.f03 os-spec-current-cyl-modes.f03 \
		os-spec-charge-cyl-modes.f03 os-spec-phasespace-cyl-modes.f03

# particle source object source files
include cyl_modes/psource/Makefile

src +=  os-spec-push-cyl-modes.f03 os-spec-cyl-modes.f03 \
		os-neutral-cyl-modes.f03 os-antenna-cyl-modes.f03 \
		os-cathode-cyl-modes.f03 \
		os-zpulse-mov-wall-cyl-modes.f03 os-zpulse-flyfoc-cyl-modes.f03 \
		os-zpulse-planewaves-cyl-modes.f03 os-zpulse-cyl-modes.f03 \
		os-particles-cyl-modes.f03 os-simulation-cyl-modes.f03

# Add local path to vpath
vpath %.f03 cyl_modes

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif

endif
