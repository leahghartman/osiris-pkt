# Overdense algorithm
### Particle stopping, splitting, and recombination

## Overview

The "overdense" algorithm is designed (1) to avoid reflux from hot particles hitting a standard boundary by stopping them and (2) to avoid the creation of an artificial wakefield due to a large, hot macroparticle by splitting hot particles into smaller ones.

## Absorber (a.k.a., damper)
A hot beam incident on a one-cell boundary (e.g., thermal, absorbing, etc.) can build up surface fields which artificially accelerate particles back into the box.  This algorithm changes the momentum of hot particles over a defined region by instantaneously "reemitting" them from a specified thermal bath.  This prevents the field from building up at the boundary.

Particle stopping works by specifying either a multiple of the local temperature or a gamma value above which particles will be considered for stopping.  A mean free path, `mfp`, is specified, which statistically determines how likely the particle is to be stopped.  In the current implementation the absorbing region is `3*mfp` long, and most particles are stopped after a distance of `2*mfp`.  When a particle is stopped, it will be "reemitted," meaning that its momentum will be instantly changed to that of a particle from a thermal bath.  The thermal velocity is either specified by the user or set as the local (within one cell) thermal velocity.

For most applications, the default values of nearly all the input parameters should work well.  The user only need specify the `n_damp`, `damper_direction`, `orientation` (either "forward" or "backward", corresponding to an absorber at the upper or lower simulation edge, respectively), `mfp` and `damp_start_time`, if desired.  The mean free path is recommended to be between 10–100  (be careful of your normalization!), depending on the application.  Factors that may require a larger `mfp` include increasing density or energy of the particle beam, along with the duration of the simulation.  Multiple `damper` sections may be specified in the input deck for a single species, since one may wish to stop particles at multiple boundaries.  Note that the `damp_species` routine is parallelized with OpenMP.

### Standard input parameters

* Following the `diag_species` section, one or multiple new sections can be specified, titled `damper`, with the following available arguments:
  * `n_damp`, integer, default = 0. Perform the damping routine every `n_damp` timesteps.
    - Note, damping every timestep could increase the advance_deposit time by \~30%, but damping should be frequent enough that a particle traveling near the speed of light has many chances (e.g., >25) to be damped over a mean free path.
  * `damper_direction`, integer, default = 1. Simulation dimension over which the damping will occur.parameters are fixed at their "max" values.
  * `orientation`, character(\*), default = "forward".
    - "forward" - places the absorber at the upper simulation boundary
    - "backward" - places the absorber at the lower simulation boundary
  * `mfp`, real, default = -1.0. Mean free path that energetic particles travel before getting stopped. The entire absorbing region will be `3*mfp` long and placed adjacent to the appropriate simulation boundary. Most particles are stopped after a distance of `2*mfp`.
  * `damp_start_time`, real, default = -1.0. Simulation time after which to begin the damping routines.

The absorber has many parameters, and if desired, a user can completely specify and fine-tune all of them.  The absorbing region is internally defined by two positions, `x_start` and `x_full`.  Stopping begins at `x_start` and continues all the way to the simulation boundary.  The user-specified parameters `vth_mult`, `gamma_co`, `reemit_temp`, and `stopping_dist` are all given two values: `*_start` and `*_full`.  The values used for these parameters for a particle located some distance into the absorbing region is linearly interpolated between the `*_start` and `*_full` values if between `x_start` and `x_full`, or set as `*_full` if past `x_full`.  Note, the `*_full` values do not have to be larger than the `*_start` values; this is just the nomenclature.

The `mfp` input parameter will automatically set `x_start`, `x_full`, `stopping_dist_start` and `stopping_dist_full` if `custom_input = .false.`.  These values can be set by the user in the input deck if `custom_input = .true.`.  All other `*_start` and `*_full` parameters may be set by the user even if `custom_input = .false.`.  All extra available parameters are listed below for users who wish to fine-tune the absorber.

### Extra input parameters

* The following arguments can also be added to a `damper` section:
  * `if_damp`, logical, default = `.true.`. Flag to turn on damping (stopping) for this species.
  * `residual`, real, default = 1.0d-7. When `custom_input = .false.` and the `mfp` is set, this parameter is used to calculate the internal `stopping_dist_full` value such that the fraction `residual` of hot particles incident on the absorber propagate all the way to the simulation boundary.
  * `custom_input`, logical, default = `.false.`. If .true., the user has to specify `x_start`, `x_full`, `stopping_dist_start` and `stopping_dist_full`.  If .false., these parameters are set internally from the user-defined `mfp` and `residual` inputs.
  * `x_start`, real, default = 0.0. The x value at which absorption begins.
  * `x_full`, real, default = 0.0. The x value at which the absorption parameters are fixed at their "max" values.
  * `local_vth_cutoff`, logical, default = `.true.`. Define the damping cutoff as a multiple of the computed local temperature (`.true.`) or as a gamma value (`.false.`).
  * `vth_mult_start`, real, default = 6.0. Cutoff multiple of the local temperature used at the start of the damper.
  * `vth_mult_full`, real, default = 6.0. Cutoff multiple of the local temperature used past `x_full`.  Linear interpolation used between `x_start` and `x_full`.
  * `gamma_start`, real, default = -1.0. Cutoff gamma value used at the start of the damper.
  * `gamma_full`, real, default = -1.0. Cutoff gamma value used past `x_full`.  Linear interpolation used between `x_start` and `x_full`.
  * `reemit_with_local_temperature`, logical, default = `.true.`. Reemit particles from a thermal distribution with temperature equal to that computed at the local cell (`.true.`) or that given by the `reemit_temp` parameters (`.false.`). Results are usually much better with this set to .true.
  * `reemit_temp_start`, real, default = -1.0. Reemission temperature of the particles used at the start of the damper.
  * `reemit_temp_full`, real, default = -1.0. Reemission temperature of the particles used past `x_full`.  Linear interpolation used between `x_start` and `x_full`.
  * `fourth_root_temperature`, logical, default = `.false.`. Use the fourth root to calculate the temperature.  This parameter affects both `local_vth_cutoff` and `reemit_with_local_temperature`.  Tests show that this can usually be left as .false.
  * `stopping_dist_start` , real, default = 0.0. The (target) stopping distance (usually 2\*mean_free_path) used at the start of the damper.
  * `stopping_dist_full`, real, default = 0.0. The (target) stopping distance (much smaller than mean_free_path) used past `x_full`.  Linear interpolation used between `x_start` and `x_full`.
  * `if_cold_perp`, logical, default = `.false.`. Reemit thermally in all three dimensions (`.false.`) or only in the physical directions of the simulation (`.true.`).
  * `linear`, logical, default = `.true.`. Whether to damp based on the velocity magnitude (`.true.`) or based on the hazard function (`.false.`). For most applications, this should be set to .true.
  * `stop_dir`, character(\*), default = "all". Setting this to "toward boundary" could help provide a more accurate return current in some cases.
    - "all" - stops energetic particles traveling in any direction
    - "toward boundary" - stops only energetic particles that are traveling in the direction of the boundary (as determined by `orientation`).
    - "away from boundary" - stops only energetic particles that are traveling away from the boundary (as determined by `orientation`).
  * `reemit_dir`, character(\*), default = "all". Setting this to "away from boundary" could potentially help the reemitted particles supply the correct return current in some cases.
    - "all" - reemits stopped particles traveling in any direction
    - "toward boundary" - only reemits stopped particles traveling in the direction of the boundary (as determined by `orientation`).
    - "away from boundary" - only reemits stopped particles traveling away from the boundary (as determined by `orientation`).
  * `use_local_vth_func`, logical, default = `.false.`. Whether to use (`.true.`) or not use (`.false.`) a function that specifies where to/not to calculate the local thermal velocity (for `local_vth_cutoff` and `reemit_with_local_temperature`).
    - The accurate calculation of the local thermal velocity relies on having a large number of particles in each cell.  If an absorber is used along a boundary where a portion of the plasma in at near-vacuum (i.e., the majority of particles are hot beam particles), the local thermal velocity may be very high, and hot particles may get through to the boundary. A user may desire thermal velocity calculation in one region of the plasma and an absolute cutoff (specified by `gamma_start`, `gamma_full`, `reemit_temp_start` and `reemit_temp_full`) in another region. If so desired, set `use_local_vth_func = .true.`.
  * `local_vth_func`, character(\*), default = "NO_FUNCTION_SUPPLIED!". Function of space used when `use_local_vth_func = .true.` to determine where to calculate the local thermal velocity (>0) or use specified values (<=0).
    - For example, if `local_vth_func = "if( x1>=0, 1, -1 )"`, then for `x1>=0` the local thermal velocity will be computed at each cell for stopping/reemitting particles, and for `x1<0` the stopping/reemission of particles will be determined by the specified energies/temperatures.

## Splitting and recombination
When a macroparticle with high energy drifts through a background plasma, it will lose energy due to the generation of a wake behind it.  The energy loss is proportional to the charge of the macroparticle, among other things (see J. May *et al.*, Phys. Plasmas, vol. 21, no. 5, p. 052703, 2014).  Particle splitting allows a user to specify how many times to split a given macroparticle above a certain energy.  When splitting a particle, the `if_push` option will impart a small fraction of the particle's transverse momentum (set by `boost`) to the two new particles, but in opposite directions.  Without enabling either `if_push` or collisions, the two new particles would occupy the exact same trajectory as the single macroparticle.  With the push, momentum is conserved, but energy is not.

A given particle will only be split once per splitting event.  It will be split into two particles, each with half the charge of the original.  The `splitting_gammas` array determines how many times particles of certain energies may be split.  For example, setting `splitting_gammas(1:2) = 1.5,2.0,` will not split particles with gamma below 1.5, will split particles once between 1.5-2.0, and split particles twice above 2.0.  The length of the specified `splitting_gammas` array determines the maximum number of times a particle may be split.  Note, if there is a density ramp, particles from the lower-density region will be treated as having been split already, as determined by their charge.

Often splitting is used in conjunction with the absorbing algorithm, which can result in a very large number of particles deposited in the stopping region.  Particle recombination allows a user to define a certain region over which the particle buffer will be searched for particles of similar charge below a certain threshold (set by `unsplit_charge_density` and `unsplit_charge_ppc`).  A recombination event merges two particles into one, preserving momentum but not energy.

To find and merge particles of like charge, the particles on each node are binned into a grid over the region specified by `recomb_box_bounds` with cell size `recomb_cell_dx`.  When a particle with a small charge is found, that particle is stored in a recombination cell.  When another particle that has been split the same number of times is found in the same cell, the particles are merged and the bin is emptied.  This is done over a loop through the particles, meaning not all possible merges will occur on a given loop.  Generally speaking, `recomb_cell_dx` should be relatively small (usually one or two simulation cells wide) to increase the number of particles recombined over one loop.

A word on the `unsplit_charge_ppc` parameter.  We desire a way to specify the charge of an unsplit particle.  However, the exact floating point representation of this charge may not be obvious (e.g., density = 1.0, num_par_x = 3).  For this case, the user would set `unsplit_charge_density = 1.0` and `unsplit_charge_ppc = 3`, and the code will compute the exact charge of an unsplit particle.  The code automatically grabs the number of particles per cell to populate the `unsplit_charge_ppc` parameter, so the user should only have to set `unsplit_charge_density`.

Due to the implementation of the splitting and recombination algorithms, these functions are not parallelized with OpenMP.  Thus the user may experience a significant slowdown if using splitting/recombination with many threads.

### Input parameters

* Following the `diag_species` section (and any `damper` sections if applicable), a new section can be specified, titled `splitting`, with the following available arguments:

  *splitting parameters*
  * `n_split`, integer, default = 0. Perform the splitting routine every `n_split` timesteps.
  * `unsplit_charge_density`, real, default = `density` obtained from `profile`. The density where particles are considered to be unsplit (should not depend on `num_par_x`).
  * `unsplit_charge_ppc`, integer, default = `product(num_par_x(1:p_x_dim))`. Number to divide `unsplit_charge_density` by to compute the actual charge of an unsplit particle. Should not need to be modified.
  * `splitting_gammas(:)`, real, default = -1.0. Array containing the gammas separating regions of no splitting, one split, two splits, etc.  For example, `splitting_gammas(1:2) = 1.5,2.0,` will not split particles with gamma below 1.5, will split particles once between 1.5-2.0, and split particles twice above 2.0.
  * `if_push`, logical, default = `.true.`. Flag to push the particles apart (slightly change transverse momentum) after splitting.  Should be set to `.true.` for nearly all applications.
  * `boost`, real, default = 0.01. Fraction of transverse momentum to add/subtract to the two particles.
  * `random_boost`, logical, default = `.true.`. Flag to multiply boost by a random number (between 0.5 and 1.5) before push.
  * `use_min_x`, logical, default = `.false.`. Flag to only perform splitting to the right or left of a certain value of `x`.
  * `min_x`, real, default = -inf. Value of `x` past which splitting will be performed.
  * `forward_of_min_x`, logical, default = `.true.`. Split in front of (`.true.`) or behind (`.false.`) `min_x`.
  * `dir_min_x`, integer, default = 1. Direction in which to check `min_x`.

  *recombination parameters*
  * `if_recombine`, logical, default = `.false.`. Flag to turn on recombination for this species. Note, recombination cannot be enabled without splitting.
  * `n_recombine`, integer, default = 25. Perform the recombination routine every `n_recombine` timesteps.  This value is often set to \~10.
   `recombine_start_time`, real, default = -1.0. Simulation time after which to begin the recombination routines.
  * `n_recombine_regions`, integer, default = 1. Number of recombination regions in use in this simulation. The maximum number is 6 (e.g., one at each edge of a 3D simulation).
  * `recomb_box_bounds(x_dim,2,n_recombine_regions)`, real, default = 0.0. Region (in simulation units) over which particles will be considered for recombination.  Second dimension corresponds to lower/upper bound, and the third to different regions.  Note, recombination regions can overlap, but this will likely not increase performance.
  * `recomb_cell_dx(1:x_dim,n_recombine_regions)`, real, default = 0.0. Size of a cell for binning particles in the recombination event.  Usually best to set near `dx` or `2*dx`.
