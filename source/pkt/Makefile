# This is meant to be included by the osiris Makefile
ifndef DISABLE_PKT
FPPF += -DENABLE_PKT

src += os-dep-density-pkt.f03 

# QED algorithm source files
src += os-pkt.f03 os-spec-define-pkt.f03 \
  os-photons-define-pkt.f03 \
	os-photons-pkt.f03 os-spec-push-pkt.f03 \
  os-particles-pkt.f03 os-simulation-pkt.f03

# Add local path to vpath
vpath %.f03 pkt

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	@$(MAKE) -C ..

clean:
	@$(MAKE) -C .. clean

endif

endif
