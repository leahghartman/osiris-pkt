# This is meant to be included by the osiris Makefile
ifndef DISABLE_QEDCYL
# Only proceed if the required sub-modules are enabled
ifndef DISABLE_CYLMODES
ifndef DISABLE_QED

FPPF += -DENABLE_QEDCYL

# qed-cyl algorithm source files
src += os-bremsstrahlung-qed-cyl.f03 os-trident-coul-qed-cyl.f03 os-spec-define-qed-cyl.f03 os-betheheitler-qed-cyl.f03 os-photons-define-qed-cyl.f03  \
		os-spec-diagnostics-qed-cyl.f03 os-photons-diagnostics-qed-cyl.f03 \
		os-spec-push-qed-cyl.f03 os-photons-qed-cyl.f03 \
		os-particles-qed-cyl.f03 os-simulation-qed-cyl.f03


# Add local path to vpath
vpath %.f03 qed-cyl

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