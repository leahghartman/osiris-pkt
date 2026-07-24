# This is meant to be included by the osiris Makefile
ifndef DISABLE_RAD
FPPF += -DENABLE_RAD

src +=  os-radio-detector.f03 os-radio-sphedet.f03 os-radio-cartdet.f03 \
		os-radio-detlist.f03 os-radio-species.f03 os-radio-particles.f03 \
		os-radio-simulation.f03

# Add local path to vpath
vpath %.f03 radio
vpath %.f90 radio

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	@$(MAKE) -C ..

clean:
	@$(MAKE) -C .. clean

endif

endif
