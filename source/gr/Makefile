# This is meant to be included by the osiris Makefile
ifndef DISABLE_GR
FPPF += -DENABLE_GR

src += os-geometry-gr.f03 os-emf-define-gr.f03 \
       os-emf-mur-gr.f03 os-emf-pec-gr.f03 \
	   os-emf-star-gr.f03 os-emf-axi-gr.f03 \
	   os-emf-boundary-gr.f03 os-emf-gridval-gr.f03 \
	   os-emf-gr.f03 \
		 os-species-define-gr.f03 os-photons-define-gr.f03 os-psource-gr.f03 \
		 os-species-charge-gr.f03 os-species-current-gr.f03 os-species-push-gr.f03 \
		 os-species-gr.f03 os-species-diagnostics-gr.f03 os-species-group-gr.f03 os-species-emit-gr.f03 \
		os-particles-define-gr.f03  os-particles-gr.f03 os-simulation-gr.f03

# Add local path to vpath
vpath %.f03 gr

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif

endif
