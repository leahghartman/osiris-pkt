# ZPULSE module

The new version of the **zpulse** module separates the different zpulse types into different objects, making it easier to extend this module and add new features while keeping the core components intact.

Defining the type of pulse to use is done through the choice of name of the zpulse section in the OSIRIS input file, e.g. a standard zpulse (instantaneous injection inside the simulation box) is done with a `zpulse{}` input file section, and a wall type zpulse (continuous injection from a simulation wall) is done through a `zpulse_wall{}` section. Here's an example of a complete input file section for point source zpulse:

```
zpulse_point
{
	a0             = 3.,               ! Normalized vector potential
	omega0         = 10.0,             ! Laser frequency in simulation units

	point_pos(1:2)  = 10.0, 0.0,       ! Position of the point source

	tenv_type       = "sin2",          ! Type of temporal envelope
	tenv_rise       = 2.82,            ! rise time
	tenv_fall       = 2.82,            ! fall time
	tenv_flat       = 0.0,             ! flat time
}
```

An OSIRIS input file can have an arbitrary number of **zpulse** sections, combining multiple types in no particular order.

## Available zpulse types

The following **zpulse** types have been implemented:

* **zpulse** - Instantaneous injection inside the simulation box
* **zpulse\_wall** - Continuous injection from a simulation wall
* **zpulse\_mov\_wall** - Continuous injection from a moving injection plane (line)
* **zpulse\_point** - Continuous injection from a point source
* **zpulse\_planewaves** - Exact injection of all field components, either instantaneous or continuous (includes flying focus implementation)
* **zpulse_flyfoc** - Flying focus laser fields using Arbitrarily Structured Laser pulses. These are based on the solutions to Maxwell's equations discussed in Pierce et al., Physical Review Research 5, 013085.

Selecting them in the input file requires setting the corresponding section to the same name.


## Example of a zpulse_flyfoc ASTRL flying focus

The below example shows how to initialize a laser pulse with a slice-dependent focal point (specied by `s0_math_func`, intensity (specified by `a0_math_func`, and local frequency (specified by `phase_derivative_math_func`). The `if_launch_t0` parameter makes the code launch the pulse fully at `t=0` instead of using the wall injection algorithm of `zpulse_mov_wall`. The envelope of this laser pulse propagates forward superluminally.

```
zpulse_flyfoc
{
	a0 = 1,
	if_launch = .true.,
	if_launch_t0 = .true.,
	omega0 = 590.4,
	pol_type = 0,
	pol = 0.0,
	if_interp_env = .false.,
	env_type = "astrl_analytic",
	s0_math_func = "26.0*(xi+5.35)",
	a0_math_func = "if( xi > (-3.1), 0, 
	                    if( xi > (-3.3), 1 * sin(1.5707 * (xi+3.1) / (0.2))^2, 
			        if( xi > (-5), 1,
				    if( xi > (-5.2), 1 * sin(1.5707 * (xi+5.2) / (0.2))^2, 0 )
				  )
		              )
			  )",
	w0_math_func = "0.16",
	xperp0_math_func = "0", 
	pol_math_func = "0",
	if_apply_divcorr = .false.,
	phase_math_func = "0.5 * 15*(xi+4.1)^2",
	phase_derivative_math_func = "15*(xi+4.1)",
}
```

## Changes to the OSIRIS input file

This new version is not fully compatible with previous input files because:

* There no longer exists a `type` parameter, since the type of **zpulse** is now defined through the section name.
* Each **zpulse** type will only accept the parameters for that particular type, e.g., point source types will not accept `per_envelope` parameters.
* **zpulse** types using continuous injection now use `tenv_*` parameters to specify the temporal envelope profile (rather then `lon_*` parameters)

Existing input files with type other than the default type will most likely not work. However, upgrading existing input files is straightforward, requiring only minimal changes. Here's an example of a wall type zpulse (zpulse_wall):

```
zpulse_wall
{
	a0             = 3.,
	omega0         = 10.0,

	tenv_type       = "sin2",
	tenv_rise       = 2.82,
	tenv_fall       = 2.82,
	tenv_flat       = 0.0,

	per_type       = "gaussian",
	per_w0         = 2.82,
	per_focus      = 20.0,
}
```

Only 3 changes were required:

* Renaming the section from `zpulse` to `zpulse_wall`
* Removing the `type = "wall"` parameter
* Renaming the `lon_*` sections as `tenv_*`
