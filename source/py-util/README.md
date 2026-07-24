# Python–Fortran interface

## Overview

The files in this `source/py-util` directory enable you to call any custom Python function as defined in a file.  These functions can take in data directly from Fortran (a copy is made), and the Python results can then be copied back directly into Fortran.

In order to simplify the Fortran/C code necessary for this implementation, all Python functions must only have one argument, which is a dictionary that we call `STATE`.  The elements of `STATE` are set in Fortran before the function call, Python performs some operation on `STATE`, and then the modified `STATE` variables can be read back after the function call.  The keys to the `STATE` dictionary must be consistent between both the Fortran and Python calls.

## Required for compilation and running

Make sure to read over this entire section.

To properly build and link the Python interface, something like the following should be defined in the configure file:

```
CONDA_ROOT = /path/to/anaconda3
PY_FCOMPILEFLAGS = -I$(CONDA_ROOT)/include/python3.9
PY_FLINKFLAGS = -Wl,-rpath,$(CONDA_ROOT)/lib -L$(CONDA_ROOT)/lib -lpython3.9
# Or maybe these linking flags
# PY_FLINKFLAGS = -Wl,-rpath,$(CONDA_ROOT)/lib $(CONDA_ROOT)/lib/libpython3.9.dylib
```

The `PY_FCOMPILEFLAGS` flag is what signals to the compiler that the Python interface is desired.  The version number of the `python3.9` folder and library should be updated as needed.

The Python calls may crash if the `-ffpe-trap=invalid,zero,overflow` flag is enabled (usually in a debug profile), so this flag may need to be removed or set to `-ffpe-trap=none`.

Finally, when the OSIRIS executable is run, the `PYTHONPATH` variable should be set to the directory (directories) containing the custom Python script(s) that will be called from OSIRIS.  More on how these scripts are called below.

## Python initialization routines

Most of you reading this are probably less interested in writing your own Python interfaces in OSIRIS and instead interested in using some of the Python initialization routines that we have built in.  We will describe those routines here, but skip down to the section "Calling Python functions" to learn more about writing your own code.

The current Python initialization routines include
* Initializing the E and B fields
* Setting external E and B fields
* Setting the exact momentum of each particle
* Setting the thermal velocity as a function of space (which is afterward sampled from a Gaussian)
* Setting the fluid velocity as a function of space
* Setting the species density as a function of space

Each of these use cases utilizes new input deck parameters and has specific keys to the `STATE` dictionary.  These are explained with examples below.  Note that in all cases, input deck parameters with names ending in `py_mod` and `py_func` have been added.  The `py_mod` parameter is the name of the Python script to be called, e.g., `py_mod = "py-script"`.  Make sure the folder containing the named script is in the `PYTHONPATH`.  The `py_func` parameter is the name of the function within the Python script to be executed, e.g., `py_func = "func_in_my_file`.  In this way, multiple functions can be defined in the same script.  Each MPI rank will make a call to the given Python functions.

Note that all arrays passed through the `STATE` dictionary will have their indices reversed in Python compared to Fortran.

A working input deck that performs most of these examples can be found in `decks/py-util/os-stdin-py-util`, along with the Python script `decks/py-util/py-script.py`.

### Setting initial and external E and B fields

To initialize fields via Python, the desired elements of `type_init_b` and `type_init_e` of the `el_mag_fld` section should be set to `"python"`.  The parameters `init_py_mod` and `init_py_func` then refer to the script and function to be called.
* In addition, the parameter `init_move_window` has been added regarding Python-initialized fields.  Calls to Python are expensive, and in a moving-window simulation, the field values are initialized every time step that the simulation domain moves.  If the initialized fields are meant to only be performed at the very first time step and not at the front edge of the moving window, then the user should set `init_move_window = .false.`.  This will avoid calling Python each time the simulation domain moves.

Similarly, to set external fields via Python, the desired elements of `type_ext_b` and `type_ext_e` of the `el_mag_fld` section should be set to `"python"`.  The parameters `ext_py_mod` and `ext_py_func` then refer to the script and function to be called for the external fields.

For example, the input deck could read like this:

```
el_mag_fld
{
  type_init_b(1:3) = "uniform", "python", "uniform",
  type_init_e(1:3) = "python", "uniform", "uniform",
  init_py_mod = "py-script", ! Name of Python file
  init_py_func = "set_fld", ! Name of function in the Python file to call (same for all components)
  ! init_move_window = .false., ! May want to declare this for a moving-window simulation

  ext_fld = "dynamic",
  type_ext_b(1:3) = "none", "python", "none",
  type_ext_e(1:3) = "python", "none", "none",
  ext_py_mod = "py-script", ! Name of Python file
  ext_py_func = "set_fld_ext", ! Name of function in the Python file to call (same for all components)
}
```

The `init_py_func` and `ext_py_func` functions are called once for each field component of each field.  The `STATE` dictionary will be prepared with the following keys and data:
* `"data"` - The full 1D, 2D, or 3D vdf of a single E or B component, including guard cells.  **This quantity should be modified with the desired field data.**
* `"fld"` - A two-character string representing the field and component, e.g., `"e1"`.
* `"x_bnd"` - A real array of size `(2, p_x_dim)` containing the real-space coordinates of the first and last elements of the vdf `"data"` in each dimension (for the local node).  Field offset in the cell is taken into account.
* `"g_xmin"` - A real array of size `(p_x_dim)` containing the global minimum bounds of the simulation.
* `"g_xmax"` - A real array of size `(p_x_dim)` containing the global maximum bounds of the simulation.
* `"nx_p_min"` - An integer array of size `(p_x_dim)` that contains the values `grid%my_nx(1,:)`, or the integer index of the first cell of the local node in the global grid along each dimension.
* `"t"` - A real array of size 1 containing the current value of time in the simulation.

Multi-dimensional Python meshgrid arrays can easily be created from the above parameters by doing the following:

```python
  # Positional boundary data
  x_bnd = STATE["x_bnd"]

  # Create x arrays that indicate the position (remember indexing order is reversed)
  nx = STATE["data"].shape
  x1 = np.linspace( x_bnd[0,0], x_bnd[0,1], nx[1], endpoint=True )
  x2 = np.linspace( x_bnd[1,0], x_bnd[1,1], nx[0], endpoint=True )
  X1, X2 = np.meshgrid( x1, x2, indexing='xy' ) # Matches Fortran array indexing
```

The field values could then be set to a function of `X1` and `X2`, e.g., `STATE["data"] = 2*(X1-1) - (X2-2)`.

#### Interpolation objects

Another common way to set the fields (or other quantities) is to interpolate from some pre-defined data onto the OSIRIS grid.  For example, let's say we add a `main` function to our Python script `py-script.py` that is meant to be run once before an OSIRIS run to generate and save an interpolation object:

```python
import numpy as np
from scipy import interpolate

if __name__ == '__main__':
    # Make position arrays
    x = np.linspace(0, 10, 100)
    y = np.linspace(0, 10, 200)

    # Compute complicated functions on a grid
    # Here we keep the same indexing as the Fortran-communicated arrays, for no real reason
    val1 = 0.2 * x[np.newaxis, :] + y[:, np.newaxis] # First function
    val2 = x[np.newaxis, :] - y[:, np.newaxis] # Second function

    # Make an interpolation object for these functions (here we set the value to 0 out of bounds)
    interp1 = interpolate.RegularGridInterpolator( (y, x), val1, method='cubic', bounds_error=False, fill_value=0 )
    interp2 = interpolate.RegularGridInterpolator( (y, x), val2, method='cubic', bounds_error=False, fill_value=0 )
    np.save( "interp", np.array([interp1, interp2]) )
```

Running `python py-script.py` will then save `interp.npy` to disk.  We can then utilize this interpolation object in an OSIRIS run to set the initial fields using `init_py_mod = "set_fld"` in this manner:

```python
def set_fld( STATE ):
    # Name of the field component
    fld = STATE["fld"]

    # Positional boundary data
    x_bnd = STATE["x_bnd"]

    # Create x arrays that indicate the position (remember indexing order is reversed)
    nx = STATE["data"].shape
    x1 = np.linspace( x_bnd[0,0], x_bnd[0,1], nx[1], endpoint=True )
    x2 = np.linspace( x_bnd[1,0], x_bnd[1,1], nx[0], endpoint=True )
    X1, X2 = np.meshgrid( x1, x2, indexing='xy' ) # Matches Fortran array indexing

    # Load interpolation variables saved earlier
    interp1, interp2 = np.load("interp.npy", allow_pickle=True)

    # Use a different interpolation object based on the field component
    if fld=="e1":
        STATE["data"] = STATE["data"] + interp1((X2, X1))
    elif fld=="b2":
        STATE["data"] = STATE["data"] + interp2((X2, X1))
```

This is a powerful way to import experimental or numerical data into OSIRIS simulations.

#### Setting fields in quasi-3D geometry

When initializing or setting external fields in the quasi-3D geometry, the following additional keys will be passed within the `STATE` dictionary:
* `"mode"` - A string representing azimuthal mode number of the current field component, e.g., `"1"`.
* `"cmplx"` - A string reading either `"re"` or `"im"`, representing if the current component is the real or imaginary part.

### Setting the particle momenta

The momentum components of all particles can be specified exactly by setting `uth_type = "python"` in the `udist` section of the input deck for a given species.  The parameters `uth_py_mod` and `uth_py_func` should also be set as the corresponding Python script and function to be called, respectively.  For example, the input deck could read like this:

```
udist
{
  uth_type = "python",
  uth_py_mod = "py-script", ! Name of Python file
  uth_py_func = "set_uth", ! Name of function in the Python file to call
}
```

To instead specify the thermal velocity of each particle as a function of space, set `use_spatial_uth = .true.` and also specify `uth_py_mod` and `uth_py_func` in the `udist` section of the input deck.  For example, the input deck could read like this:

```
udist
{
  use_spatial_uth = .true.,
  uth_py_mod = "py-script", ! Name of Python file
  uth_py_func = "set_uth", ! Name of function in the Python file to call
}
```

To specify the fluid velocity of each particle as a function of space, set the `ufl_py_mod` and `ufl_py_func` parameters in the `udist` section of the input deck.  For example, the input deck could read like this:

```
udist
{
  ufl_py_mod = "py-script", ! Name of Python file
  ufl_py_func = "set_ufl", ! Name of function in the Python file to call
}
```

The `uth_type = "python"` and `use_spatial_uth = .true.` commands cannot be combined when setting with Python, because either the exact fields should be set or the thermal velocity should be set.  However, the fluid velocity can be included in either case.  See the input deck in `decks/py-util/os-stdin-py-util` as an example.

In each of the above cases, the `STATE` dictionary will be prepared with the following key:
* `"x"` - A real array of size `(p_x_dim, npart)` containing the positions of the particles.

The desired momentum array can then be created and set based on the positions `"x"`.  This array should be passed to the `STATE` array with the following key:
* `"u"` - A real array of size `(3, npart)` containing either the thermal or fluid momenta of the particles.  **This quantity should be set to the desired momentum data.**

Note that the Python functions for setting momenta will be called repeatedly, as the particle initialization routines process `p_cache_size` particles at a time.  See the functions `set_uth` and `set_ufl` in `decks/py-util/py-script.py`.

### Setting the species density

A new profile type has been added to import a uniformly gridded density profile that is then linearly interpolated at each particle position to define the density.  This is done by setting `init_type = "python"` in the `species` section of the input deck and by defining the `py_mod` and `py_func` parameters in the `profile` section of the input deck.  For example, the input deck could read like this:

```
species
{
  num_par_x(1:2) = 2, 2,
  init_type = "python",
}

udist
{
}

profile
{
  py_mod = "py-script", ! Name of Python file
  py_func = "set_density", ! Name of function in the Python file to call
}
```

In this case, the `STATE` dictionary does not provide any initial data to the Python function.  Rather, the `STATE` dictionary **should be set** with the following keys:
* `"nx"` - An integer array of size `(p_x_dim)` indicating the size of the density grid in each dimension.
* `"xmin"` - A real array of size `(p_x_dim)` indicating the lower `x` position of the density grid in each dimension.
* `"xmax"` - A real array of size `(p_x_dim)` indicating the upper `x` position of the density grid in each dimension.
* `"data"` - A regularly spaced real array of the same size as specified by `"nx"`.  **This quantity should be modified with the desired density data.**

Any positions outside of the region marked by `"data"` are given the same density as the outermost point that is closest.  See the `set_density` function in `py-script.py` for an example on how to do this.

## Calling Python functions

Suppose that you have a Python script titled `my_script.py` with a function `add_const` that takes one dictionary as an argument.  First off, OSIRIS does not need to know about this file at compile time; you only need to specify the path to `my_script.py` in the `PYTHONPATH` in the console at runtime, as mentioned above.

Let's say that in `add_const`, you want to read in a 1D double array, `x`, from Fortran and add a double, `y`, as supplied by Fortran to each element in the array.  Congratulations for wanting to do this in Python; it should be much more efficient than doing it in Fortran.  Your Python file, `my_script.py`, should look something like the following:

```python
# Outside of the function definition, any included code only runs once at the
# initialization. You could, e.g., load some external stuff here

def add_const( STATE ):
  x = STATE["x"]
  y = STATE["y"]
  x = x + y
  STATE["x"] = x
```

You can see that we unpack the `STATE` variables, edit one of them, and then store it back in the `STATE` dictionary.  The final version of `STATE` will then be accessible to Fortran after the function call.

Now, you'll need to call the Python interface from Fortran.  You could write a Fortran subroutine that looks something like the following:

```
subroutine add_const( x, y )

  use m_parameters
  use m_callpy ! This is they key module for calling Python

  implicit none

  real(p_double), dimension(:), intent(inout) :: x
  real(p_double), intent(in) :: y

  ! Set the STATE dictionary with the appropriate keys/values (copies are made)
  call set_state( "x", x )
  call set_state( "y", (/y/) )

  ! Call your python function as defined in your python file
  call call_function( "my_script", "add_const" )

  ! You can also call built-in functions of STATE, e.g., printing the entire dictionary
  call call_function( "builtins", "print" )

  ! Get any STATE dictionary values back
  call get_state( "x", x )

end subroutine add_const
```

Make sure that the first argument to the `call_function` function is the name of your Python script; the second argument is the name of the function you wish to call.  Also, the name of variables do not need to be the same between Fortran and Python, i.e., it's perfectly fine to write `call set_state( "arr", x )` (as long as you then reference `STATE["arr"]` in Python).

For a working example of this kind, see the `os-callpy-example.f03` and `py-script.py` files in the `source/py-util` directory.

### Some things to be aware of

Currently, there is no interface to add a lone `double` to the `STATE` variable, so we make it an array of size 1 in Fortran.  Make sure you remember this when you get to Python land (the above example actually works as intended).

The `STATE` variable persists for the duration of the OSIRIS execution.  So consider writing over dictionary keys with large arrays if you end up running out of memory somehow.

Speaking of memory, any call to `set_state` from Fortran immediately makes a copy of whatever variable you specify.  Overhead can be reduced by limiting the number and size of arrays that are passed to Python via `set_state`.

Data in Fortran is stored in the column-major format, whereas the row-major format is the default in numpy.  This means that your indexing will be reversed in numpy compared to Fortran.

## Summary

In this way, you can essentially define a Python function with any number or type of arguments, but the interface to each function is the same (requires much less messy Fortran/C code).  Happy coding!!

## Licensing and Citation

Everything in the `source/py-util` directory is licensed under the OSIRIS GNU Affero General Public License except for the contents of the subdirectory `source/py-util/call_py_fort`.  The `callpy.f03` file was adapted from the `src/callpy_mod.f90` file from the [call_py_fort](https://github.com/nbren12/call_py_fort) repo by nbren12, which is distributed under the Apache License 2.0 (see license placed in `source/py-util/call_py_fort` and `codegen/py-util`) and is referenced at the following Zenodo DOI: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7779572.svg)](https://doi.org/10.5281/zenodo.7779572).  The `plugin.c` and `plugin.h` files were generating from running `python builder.py` in the `codegen/py-util` folder.  Everything in the `codegen/py-util` directory is also distributed under the Apache License 2.0.
