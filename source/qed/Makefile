# This is meant to be included by the osiris Makefile
ifndef DISABLE_QED
FPPF += -DENABLE_QED

# Particle merging
src += os-mergedel.f03

# QED algorithm source files
src += os-qed.f03 os-qed-coulomb.f03 os-bremsstrahlung.f03 os-trident-coul.f03 os-spec-define-qed.f03 \
	os-betheheitler.f03 os-photons-define.f03 os-spec-diagnostics-qed.f03 \
	os-photons-diagnostics.f03 os-photons.f03 os-spec-push-qed.f03 \
	os-qed-collide.f03 os-particles-qed.f03 os-simulation-qed.f03


# Add local path to vpath
vpath %.f03 qed

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	@$(MAKE) -C ..

clean:
	@$(MAKE) -C .. clean

endif

endif