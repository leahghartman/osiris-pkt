# This file has been very slightly modified from the src/builder.py file in the
# call_py_fort repo by nbren12.

# file plugin_build.py
import os
import cffi

ffibuilder = cffi.FFI()

header = """
extern int set_state_py(char *, double *, int*, int*, int*);
extern int set_state_char(char *, char *);
extern int get_state_char(char *, char *, int *);
extern int get_state_py(char *, double *, int*);
extern int set_state_1d(char *, double *, int*);
extern int set_state_scalar(char *, double *);
extern int call_function(char *, char *);
"""

with open("plugin.h", "w") as f:
    f.write(header)

ffibuilder.embedding_api(header)

ffibuilder.set_source(
    "my_plugin",
    r"""
    #include "plugin.h"
""",
)

# Default to opening callpy.py instead of taking in an argument
with open("callpy.py") as f:
    ffibuilder.embedding_init_code(f.read())

ffibuilder.emit_c_code("plugin.c")
# ffibuilder.compile(target="libplugin.so", verbose=True)
