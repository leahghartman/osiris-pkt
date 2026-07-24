# This is meant to be included by the osiris/py-util Makefile

# Python utility functions
src += callpy.f03

ifdef PY_FCOMPILEFLAGS
src += plugin.c
endif

# Add local path to vpath
vpath %.f03 py-util/call_py_fort
vpath %.c py-util/call_py_fort

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
  @$(MAKE) -C ..

clean:
  @$(MAKE) -C .. clean

endif