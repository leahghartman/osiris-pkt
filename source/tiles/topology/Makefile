# Topology source files
src += os-terminal-chars.f03
src += os-hilbert.f03
src += os-snake.f03

# Add local path to vpath
vpath %.f03 tiles/topology

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif
