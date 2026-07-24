## Additional input parameters
* The following parameters have been added to `udist`, when using `uth_type = "relmax boosted"` (relativistic Maxwell-Juettner distibution with a boost):
  * `use_spatial_relmax_boost`, logical, default = `.false.`.
  If `.true.`, then the relativistic boost and relativistic temperature can have a spatial dependance, which are specified by `spatial_beta_boost` and `spatial_relmax_T`, respectively.
  If `.false.` then both the relativistic boost and relativistic temperature are constant in space and are specified by `relmax_beta_boost` and `relmax_T`.
  * `spatial_relmax_T`, character(*), default = `0`.
  This is a mathematical function defining the temperature of the relativistic Maxwell-Juettner distribution (which may vary in space). Importantly, this temperature is defined in the rest frame of the particle species (i.e., where the net momentum/motion of this species in any direction is zero). The temperature in the rest frame is related to the temperature in the simulation frame (or lab frame) by $T_{rest} = \gamma(\beta_{boost}) T_{sim}$.
  * `spatial_beta_boost`, character(*), default = `0`.
  This is a mathematical function defining beta ($v_{boost}/c$) of the distribution in the simulation frame (which may vary in space).
  * `relmax_boost_dir`, integer, default = `-1`.
  This sets the direction of the boosted velocity (only relevant for beta_boost). In the current implementation, the boost can only be applied along one of the principal axes of the simulation (1,2 or 3).
  * `relmax_beta_boost`, real, default = `0.0`.
  Specifies constant boost velocity $\beta  = v_{boost}/c$ of the distribution function.
  * `relmax_T`, real, default = `0.0`.
  Specifies constant relativistic temperature of the distribution. Note that this temperature is defined in the rest frame of the plasma.

  Below is an example of the initialization of a relativistic Maxwell-Juettner distribution with a spatially varying temperature and beta_boost:
    ```
    udist
    {
      uth_type = "relmax boosted",
      use_spatial_relmax_boost = .true.,
      spatial_beta_boost = "0.5*exp(1-sqrt(x1^2 + x2^2))",
      spatial_relmax_T = "(1. / sqrt(x1^2 + x2^2)) / sqrt( 1. - ( 0.5*exp(1-sqrt(x1^2 + x2^2)) )^2 )",
      relmax_boost_dir = 3,
    }
    ```

* A new species pusher called `"gravity"` has been implemented that adds a uniform acceleration to the dudt of each particle in order to mimic the presence of a graviational field. To use this pusher, set `push_type = "gravity"` for each species section. The parameters of this pusher are as follows:
  * `gravity`, real(3), default = `0.0`.
  This parameter sets the uniform gravitational acceleration experienced by a given species. The components of the graviational field along the principal cartesian axes need to be specified. The graviational acceleration is expressed in units of $c\omega_{p}$.

  Below is an example of a species section that uses the gravity pusher:
    ```
    species
    {
      name = "electrons",
      push_type = "gravity",
      num_par_max = 1000000,
      rqm=-1.0,
      num_par_x(1:2) = 20, 20,
      add_tag = .false.,
      free_stream = .false.,
      gravity(1:3)= -0.005, 0.0, 0.0,
    }
    ```

* The following parameters have been added to the `diag_species` section under `rep_udist` and specify grid-resolved fluid and pressure tensor diagnostics, including spatial/time averaging, lineouts, etc., as described in the grid diagnostics section.
  * `vfl1`, `vfl2`, `vfl3` - Fluid velocities, defined as  $v_{{fl}_i} \equiv \frac{1}{n_s}\int {d^3u f_s(x,u) u_i/ \gamma}$, where $f_s(x,u)$ is the distribution function of species s, and i is the ith direction.
  * `P11`, `P22`, `P33` - Diagonal components of the momentum flux tensor, $P_{ii} \equiv \int {d^3u f_s(x,u) u_iu_i/ \gamma}$.
  * `P12`, `P13`, `P23` - Off-diagonal commponents of the momentum flux tensor, $P_{ij} \equiv \int {d^3u f_s(x,u) u_iu_j/ \gamma}$.
  * `nufl1`, `nufl2`, `nufl3` - Fluid momentum density, defined as $n u_{{fl}_i} \equiv \int {d^3u f_s(x,u) u_i}$.