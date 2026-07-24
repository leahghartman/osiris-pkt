# this is ugly but using I don't understand ($or..) in make
ifndef DISABLE_PGC
        USE_TRIDIAG = "yes"
endif
ifndef DISABLE_SHEAR
        USE_TRIDIAG = "yes"
endif


ifdef USE_TRIDIAG
        src += os-tridiag.f03
endif

ifndef DISABLE_PGC
FPPF += -DENABLE_PGC

# PGC algorithm source files
src += os-emf-pgc-define.f03 os-emf-diag-pgc.f03 \
                         os-pgc-laser.f03 os-emf-pgc.f03 \
                         os-pgc-interpolate.f03 os-pgc-deposite.f03 \
                         os-spec-pgc.f03 os-neutral-pgc.f03 os-particles-pgc.f03 \
                         os-simulation-pgc.f03

endif

# Add local path to vpath
vpath %.f03 pgc

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
        $(MAKE) -C ..

clean:
        $(MAKE) -C .. clean

endif
