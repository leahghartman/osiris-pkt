# This is meant to be included by the osiris Makefile

ifneq ("$(ENABLE_PY_UTIL)","")
  $(error ERROR Specifying ENABLE_PY_UTIL in configuration file is now deprecated. \
    Instead, specify only PY_FCOMPILEFLAGS and PY_FLINKFLAGS.)
endif

# call_py_fort
include py-util/call_py_fort/Makefile

ifdef PY_FCOMPILEFLAGS
src += os-callpy-example.f03
endif

# Add local path to vpath
vpath %.f03 py-util

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
  @$(MAKE) -C ..

clean:
  @$(MAKE) -C .. clean

endif