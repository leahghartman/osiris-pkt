# This is meant to be included by the osiris Makefile
ifndef DISABLE_TILES
FPPF += -DENABLE_TILES

# Tile topology files
include tiles/topology/Makefile

# Tiles algorithm source files
src += os-tconf.f03 os-nconf-tiles.f03 \
		os-grid-tiles.f03 os-vdf-comm-tiles.f90 \
		os-wall-comm-tiles.f90 os-spec-comm-tiles.f90 os-spec-rawdiag-tiles.f90 \
		os-spec-phasespace-tiles.f90 os-spec-loadbalance-tiles.f90 \
		os-bnd-tiles.f03 os-particles-tiles.f90 os-simulation-tiles.f03 \
		os-simulation-group.f03 os-group-comm.f03 \
		os-group-boundary.f03 os-dutil-tiles.f03 os-dlb-group.f03

# Add local path to vpath
vpath %.f90 tiles
vpath %.f03 tiles

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif

endif
