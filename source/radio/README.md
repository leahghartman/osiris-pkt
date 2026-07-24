# RaDiO:wiki

In this page you will find everything there is to know about running and developing the radiation tools for OSIRIS.

# User Guide

## Input File

The input file should follow the same format as a standard OSIRIS input deck,
however, there are some minor changes  to take into account:

### Simulation mode
The first thing to specify in the input file is the simulation mode:
``` 
simulation
{
  algorithm="radiat",
}
```

### Number of radiative species
Radiative species are the ones whose radiation is tracked. It is possible to run
simulations with a mixture of radiative and non-radiative species.

The number of radiative species is defined by `num_radiat` and the total number 
of species is given by `num_species+num_radiat`

``` 
particles
{
  num_species=1,
  num_radiat=1,
}
```

In this example the simulation will have 2 different species: a regular one
and a radiative one. 

It is also important to note that the radiative species come first in the
input deck.

#### Neutral species

Recent changes made it possible  to track radiation from electrons from ionized neutrals
through the `num_neutral_rad` and `num_neutral_mov_ions_rad` parameters

``` 
particles
{
  num_species=1,
  num_neutral=1,
  num_neutral_rad=1,
  num_neutral_mov_ions_rad=1,
}
```

In this case, there are 4 species: a regular one and three neutrals, one of which is radiative and one of which is radiative with moving ions.
Please note that all radiative species come before their regular counterparts in the input deck.

### Detector

One or more detector may be associated with each radiative species. The detector
is region of space and time where the radiation from that species is tracked and
it is the last thing to be specified in the input deck. 

A detector may either be `spherical` or `cartesian` as explained in the "Implementation"
section. 
Defining the type of detector to use is done through the choice of name of the 
detector section in the OSIRIS input file, e.g. a spherical detector is done with a `spherical{}`
input file section, and a cartesian type detector is done through a `cartesian{}` section. 
It is also possible to perform a selection of the radiative particles through a `selection{}`
section, with similar syntax to the `raw` selection rules and also with the possibility to use tags.
However, note that the `rad_math_expr` variable will override and cause to ignore the `rad_gamma_limit` parameter.
Here's an example of a complete input file section for a spherical detector:

``` 
species1{...
}

profile{...
}

diag_species{...
}

!optional section
selection
{
  rad_fraction=0.01,
  rad_math_expr="x1-t>2"
  !rad_gamma_limit=2,
  !rad_tags=.True.,
  !rad_file_tags="species 1_.tags",
}

spherical{
  nf(:)   = 4096,1,1024,1,
  dmin(:) = -20,1.5,-0.1,10000,
  dmax(:) = -10,1.6,0.1,100001,
  approx  = .False.,
  comps(:)= .False.,.True.,.False.,.False.,.False.,.False.,
  ndump_fac_det=900,
}
```
An OSIRIS input file can have an arbitrary number of detector sections, 
combining multiple types in no particular order, provided that they come at the end of
the species section.

Each of the detector parameters are explained below:

#### Boundaries
The detector must be specified as a 4-Dimensional range. In order to do that, the lower and upper boundaries must be defined.
The boundaries can be specified with the following lines:
``` fortran
!lower and upper bounds of the detector: (t,x1,x2,x3)
dmin(:)=tmin, x1min, x2min, x3min,
dmax(:)=tmax, x1max, x2max, x3max,
```
In the spherical version x1, x2 and x3 correspond to $ \theta,~\phi,~r $, respectively.

#### Number of Cells
The number of cells in direction is defined with the <code>nf</code> parameter:
``` fortran
!number of cells in each direction
nf(:)=1024,512,512,1 
``` 
In case there is only one cell in a given direction, that cell will be located at the corresponding <code>dmin</code> parameter.


#### Output components
In order to select which components of the field are to be calculated. The selected components can be toggled using the <code>.True.</code> and <code>.False.</code> flags. However, in order for any B-field components to be reported, you must compile the code with the `__HAS_RAD_BFLD__` flag defined (can be turned on in `source/os-config.h`). Otherwise they will automatically be set to false.
``` fortran
comps(:)=.False.,.True.,.False.,.False.,.False.,.False.,
! (Ex,Ey,Ez,Bx,By,Bz)
``` 

#### Far-Field Approximation

The Far-Field approximation section can be toggled using the <code>approx</code> parameter. 
``` fortran
! Electric field calculation option
! .false. => Exact calculation, .true. => Far-field approximation
approx = .true.,
``` 
This approximation requires less calculations than the exact calculation so it will speed up the code when selected.

## Performance
This section contains information about how several aspects of the code affect performance.


#### Time Resolution
In each time step, for each particle the algorithm circles through all detector cells, and finds at which time slot it should be deposited. Therefore, the number of operations does not depend on the time resolution, and the only effect of increasing time resolution is the increase of the detector's size which may affect the communications' performance. However, communications are not abundant in both version of the code:
*In the post-processing version, communications only happen at the end, when all processing units share their local detectors before outputting the global detector in a single file
*In the run-time version, these same communication occur whenever the detector diagnostic is used.
#### Spatial Resolution
The number of operations per time step is given by <math>N_cN_p</math>, where <math>N_c</math> is the number o detector cells and <math>N_p</math> is the number of radiative particles. This way, increasing the number of spatial cells in the detector increases the number of operations.

