# GR Module 
 
## 1. Introduction 
OSIRIS-GR is an extension of OSIRIS developed to model the global electrodynamics of extreme environments surrounding compact objects such as neutron stars and black holes. It has been developed by Rui Torres and Fábio Cruz (IST). 

## 2. Physics 
### Field equations 
OSIRIS-GR solves a set of Maxwell’s equations modified by the presence of strong gravitational fields [1,2] driven by massive rotating compact objects,
$$\partial \mathbf{B} / \partial t = - c \nabla \times ( \alpha \mathbf{E} + \boldsymbol{\beta} \times \mathbf{B} ) \ ,$$
$$\partial \mathbf{E} / \partial t = c \nabla \times ( \alpha \mathbf{B} - \boldsymbol{\beta} \times \mathbf{E} ) - 4\pi \mathbf{J} \ ,$$
where $$\alpha$$ is the lapse function, that measures the ratio of the ticking rate of the local observer to the ticking ratio of the universal time and $$\boldsymbol{\beta}$$ is the shift vector, the velocity of the spacetime shearing in relation to the grid. Without a gravitational field, $$\alpha = 1$$ and $$\boldsymbol{\beta} = 0$$ and we recover the ordinary Maxwell’s equations. The current density in Ampère's law is also related to $$\alpha$$ and $$\boldsymbol{\beta}$$, $$\mathbf{J} = \alpha \mathbf{j} - \rho \boldsymbol{\beta}$$, where $$\rho$$ and $$\mathbf{j}$$ are the ordinary plasma charge and current densities. 

### Particle equations 
Particles’ equations of motion are also modified in strong gravitational fields. In their most general form [1,2], they become 
$$\frac{\partial \boldsymbol x}{\partial t} = \frac{\alpha}{\Gamma}\boldsymbol u - \boldsymbol \beta$$
$$\frac{\partial \boldsymbol u}{\partial t} = \frac{\alpha q}{m}\left(\boldsymbol E + \frac{\boldsymbol u}{\Gamma}\times\boldsymbol B\right) + \alpha \Gamma \boldsymbol g + \alpha \overset{\leftrightarrow}{\boldsymbol H} \cdot\boldsymbol u + \frac{\alpha}{m}\boldsymbol{\mathcal{F}}_\text{RR} - \boldsymbol{\mathcal{C}}_\text{Coord},$$
where $$\boldsymbol g$$ is the gravitational acceleration, $$\overset{\leftrightarrow}{\boldsymbol H}$$ is the gravitomagnetic tensor, $$\Gamma = \sqrt{\varepsilon +\boldsymbol u\cdot\boldsymbol u }$$ is the particle's special relativistic Lorentz factor, and $$\varepsilon = 0, 1$$ for massless and massive particles, respectively. The last two forces are, respectively, the Landau-Lifshitz radiation reaction force [3] and a non-inertial, coordinate dependent, force given by
$$\boldsymbol{\mathcal{F}}_\text{RR} = \frac{2e^3}{3mc^3}\left[\frac{e}{mc}\left(\boldsymbol E\times \boldsymbol B + \boldsymbol B \times \left(\boldsymbol B \times \frac{\boldsymbol u}{\Gamma} \right)+\boldsymbol E\left(\frac{\boldsymbol u}{\Gamma}\cdot \boldsymbol E\right)\right) - \frac{e\Gamma}{mc}\boldsymbol u \left(\left(\boldsymbol E +\frac{\boldsymbol u}{\Gamma} \times \boldsymbol B\right)^2 - \left(\boldsymbol E \cdot\frac{\boldsymbol u}{\Gamma} \right)^2\right) \right]$$
$$\boldsymbol{\mathcal{C}}_\text{Coord} = \mathcal{C}_\text{Coord}^{\hat i} = \frac{h_i}{h_k} {}^3\Gamma^i_{jk} \frac{dx^j}{dt}u^{\hat k} - \frac{h_i}{h_i}{}^3\overset{\text{NS}}{\Gamma^i_{ij}}\frac{dx^j}{dt}u^{\hat i} \ ,$$
where $${}^3\Gamma^i_{jk}$$ are the purely spatial Christoffel symbols (NS stands for not considering the sum over repeated indexes) and $$h_i$$ are the Lamé coefficients. 

## 3. Algorithms 
### Field solver 
OSIRIS-GR is currently implemented in 2D spherical coordinates ($$r, \theta$$), i.e. it assumes symmetry about the azimuthal angle $$\phi$$. The simulation grid is uniform in $$\log r$$ and $$\theta$$. The modified set of Maxwell’s equations are solved with an integral form to avoid divergences at the polar axes [e.g. 4], that transforms the curl operators as 
$$(\nabla \times \mathbf{E})_i = \left( \oint \mathbf{E} \cdot \mathrm{d} \mathcal{C}_i \right) / \mathcal{S}_i \ ,$$ 
where $$\mathcal{S}_i$$ is the area element with normal perpendicular to direction $$i$$ and $$\mathcal{C}_i$$ is a curve defining its boundary. 
This module supports multiple spacetime metrics: Minkowski, Schwarzschild and the slow rotation limit of Kerr (see section 4. for details). The properties of the adopted spacetime metric are embedded in the volume, area and line elements used to solve Maxwell’s equations given in section 2. using the integral form described above. 

### Particle pusher 
Particles are pushed with Runge-Kutta (RK) methods. Both 4th and 6th order RK methods are available and can be chosen at run time. These methods evolve the orthonormal contravariant components of the velocity and contravariant components of the position in spherical coordinates. Despite the axisymmetry assumption, all three components of the velocities and positions are advanced. To avoid singularities on the polar axes, particles are pushed in cartesian-like coordinates close to the polar axes. For simulations without GR effects, a Boris pusher is also available that advances particles’ position and velocity components in cartesian coordinates ($$x,y,z$$). 

### Current deposition 
Current is deposited using a charge-conserving scheme developed specifically for the non-uniform curvilinear grid considered in OSIRIS-GR [5]. For non-flat spacetimes, the quantity deposited is $$\mathbf{J} = \alpha \mathbf{j} - \rho \boldsymbol{\beta}$$ instead of the ordinary current density $$\mathbf{j}$$. This scheme will be detailed in a future publication.

### Boundary conditions 
The inner radial boundary is expected to mimic the compact object (either the neutron star surface or the black hole event horizon) and is thus at a finite radius $$r = r_* > 0$$. This boundary is treated as a perfect conductor for fields and can be either assumed to be static or rotating. In the latter case, the conductor is assumed to start at rest, and is linearly spun up over a given time, a condition imposed via the boundary condition. The outer radial boundary can either be assumed to open (Mur) or conducting. Both radial boundaries are open for particles. 
The simulation domain is expected to span the whole polar angle $$\theta \in [0, \pi]$$. The boundaries in this direction assume symmetry in the azimuthal direction and work effectively as reflecting boundaries for both fields and particles. 

### Strong-field effects 
Pair production mechanisms typical of extreme astrophysical settings are mimicked using phenomenological models. For the time being, only one model is available that produces new pairs whenever electrons or positrons reach a threshold Lorentz factor. Pair particles are created with an energy that is a fixed fraction of their parent particle, and with momentum parallel to their parent. 

## 4. User guide

**Important:** OSIRIS-GR should be compiled with double precision, since the calculation of sensitive functions required in this simulation mode such as $$\log(x) - \log(x+dx)$$ or $$x^p - (x + dx)^p$$ in single precision may lead to significant round-off errors.

In order to use the OSIRIS-GR module, select it as a simulation mode in the `simulation{}` section of the input file: 
``` 
simulation 
{ 
  algorithm = "gr",  !  simulation mode 
} 
``` 

Recall that the grid is spherical, which means that the first entry refers to the radial direction and the second to the poloidal direction (being $$\theta = 0$$ the north pole). In this way, we have: 
``` 
!-- spatial grid 
grid 
{ 
  nx_p(1:2)   =  1000 , 1000,  !  cell number along r and theta 
} 
!-- spatial limits 
space 
{ 
  xmin(1:2) = 1.0, 0.0,            !  minimum position along r and theta 
  xmax(1:2) = 20.0, 3.1415926536,  !  maximum position along r and theta 
} 
``` 

### Geometry
For this simulation mode, a new `geometry{}` section is required in the input file. This section allows the specification of the geometric properties of the problem. Here is an example of a complete geometry section: 
``` 
geometry 
{ 
  metric = "kerr_slow",      !  type of spacetime metric 
  rs = 0.7,                  !  normalized Schwarzschild radius 
  bprofile = "dipole",       !  type of magnetic field profile 
  bs = 3.53190d3,            !  normalized stellar polar surface magnetic field 
  omega = 0.125,             !  normalized stellar angular velocity 
  trise = 1.5,               !  normalized time for angular velocity saturation 
  classical_rr = .true.,     !  selects if using classical radiation reaction 
  b_schwinger = 8.0d4,       !  normalized Schwinger field value 
  spline_rule = "integral",  !  type of spline rule used in charge cons./current dep. 
} 
``` 

An OSIRIS input file must have only one geometry section. 
#### Available metrics 
The `metric` parameter defines the spacetime metric used in the simulation. The following metrics are available: 
- `"minkowski"` - Flat spacetime (`rs = 0`), default 
- `"schwarzschild"` - Curved spacetime around a non-rotating compact object 
- `"kerr_slow"` - Curved spacetime around a slowly rotating compact object 

If `rs = 0`, the metric is automatically chosen to be `"minkowski"` (default). The other metrics are selectable for `rs > 0.0001`. 

#### Available initial magnetic field profiles 
The parameter `bprofile` is used to define the type of magnetic field used in the electromagnetic field boundary conditions. This should be consistent with the magnetic field profile initialized in section `el_mag_fld{}`. The following magnetic field profiles are available: 
- `"monopole"` - Monopolar magnetic field configuration (only available with Minkowski metric) 
- `"dipole"` - Dipolar magnetic field configuration 

When using a curved spacetime, the dipolar field is modified to the dipolar field solution in vacuum with the corresponding selected spacetime metric. This modification includes non-dipolar field components [6]. 
The default magnetic field profile is set to `"none"`. 

#### Available spline rules 
The parameter `spline_rule` selects the rule to be used when evaluating the spline functions used to deposit the charge on the grid and, consequently, calculate the required electrical current density. The following spline rules are available: 
- `"rectangular"` - Approximate the radial integral by using the rectangular rule 
- `"integral"` - Evaluates the complete integral 

Both rules yield machine precision charge conservation in flat spacetime configurations. However, when using a curved spacetime, the integral rule is the only one that conserves charge to machine precision. For that reason, the rectangular rule must be used with care. 
The default spline rule is set to `"integral"`. When not specifying the used spline rule, the charge density diagnostic will use a simple linear spline on the mid-cell position. 

#### Other parameters 
`rs` is a parameter that selects the normalized radial position of the event horizon of the compact object. This value varies from 0 (flat spacetime - Minkowski) to 1. This value must never reach unity as this would mean that the inner boundary is indeed the event horizon. Further development is required to allow for this condition. When using the slowly rotating limit of the Kerr metric, we recommend not exceeding the 0.7 value. The default value is set to 0, in accordance with the default metric type. 

`bs` is a parameter that selects the normalized magnetic field intensity on the stellar surface and at the pole. The default value is set to 100. Recall that this should be in accordance with the initialized magnetic field amplitude on the `el_mag_fld{}` section. When using a curved spacetime metric, the field intensity is modified and `bs` is no longer the value on the stellar polar surface. The selected value on the example is such that the field intensity on the stellar polar surface is 8000 for `rs = 0.7`.

`omega` is a parameter that selects the normalized stellar angular velocity. This parameter also defines the light cylinder radius $$R_{LC} = c / \Omega$$. The radial extent of the simulation box is recommended to be at least twice the light cylinder radius. By default, this value is set to 0.2. 

`trise` is a parameter that selects the normalized time it takes for the star to reach the selected angular velocity. In this simulation mode, the star always starts with zero angular velocity and then is spun to omega via the stellar electric field boundary condition. We recommend this value to be at least 1 for a smooth initial radially outgoing electromagnetic field perturbation. The default value is set to 10. 

`classical_rr` is a boolean parameter that selects if the classical radiation reaction force is present (`.true.`) or not (`.false.`) in the simulation. This is only possible for curved spacetime metrics. By default, this value is set to `.false.`. 

`b_schwinger` is a parameter that sets the normalized value of the critical magnetic field (Schwinger magnetic field). This value should be larger than any grid value of the initialized magnetic field. We recommend it to be at least 10 times the stellar polar surface magnetic field (`bs`). The default value is set to 10^6. This parameter is only important if `classical_rr` is set to `.true.`. 

#### EM field boundary conditions
For more details on the boundary conditions, read the section on Boundary Conditions above. Here we present an example of a typical `emf_bound{}` input file section: 
```
emf_bound 
{ 
  type(1:2,1) = "star", "mur", 
  type(1:2,2) = "axisymmetric", "axisymmetric", 
} 
```

The available boundary condition types are: 
- `"pec"` or `"conducting"` - Perfect electric conductor (can be applied to the inner and/or outer radial boundaries) 
- `"star"` - Rotating magnetized conductor (can only be applied to the inner radial boundary) 
- `"axial"` or `"axisymmetric"` - Imposes the axial symmetry on both poloidal boundaries 
- `"open"` or `"mur"` - Perfect absorber for purely radial perturbations, mimicking an open radial boundary (can only be applied to the outer radial boundary) 

Recall that some of these boundary types are special, which means that they can only be adopted in only one specific boundary (see specifications above). These conditions must be specified, otherwise an error will pop up (there are no default values!). 

### Species and species groups
Apart from the typical species objects, OSIRIS-GR also supports species groups. These groups facilitate injection and/or production of copious pairs in the simulation to mimic strong-field effects. To include a new species group in a simulation, a `group{}` section should be added to the input file and followed by two species blocks (including `species{}` and the corresponding blocks for spatial and velocity distributions and boundary conditions). The `group{}` section should look like the following example:
``` 
group 
{ 
  cooling_type = "qed_red",  ! particle cooling mechanism 
  pp_gth = 40,               ! threshold Lorentz factor for pair production 
  pp_gpair = 12.8,           ! Lorentz factor of the created pair 
  pp_rmin   = 1.03,          ! lower radial boundary for pair production 
  pp_rmax   = 3.,            ! upper radial boundary for pair production 
  inj_local = "surface",     ! particle injection mechanism 
  inj_type = "fracEB",       ! charge density injection mechanism 
  ks = 0.2,                  ! fraction of the charge density to inject 
  vs = 0.0,                  ! particle velocity at injection 
  inj_limit = "fracEB",      ! injection limiting criterion 
  k_limit = 5.0,             ! injection limiting parameter 
} 
``` 

#### Available cooling types 
`cooling_type` selects the particle cooling mechanism. For the time being, only a reduced pair creation model, that produces new pairs whenever electrons or positrons reach a certain threshold Lorentz factor, is available. This reduced model can be selected by choosing `"qed_red"`. The default cooling type is set to `"none"`. 

#### `qed_red` parameters 
`pp_gth` is a parameter that defines the threshold Lorentz factor for pair production. A new pair is created when any electron or positron within the species group reaches a Lorentz factor `pp_gth`. Near the polar axes, this value is multiplied by an increasing function that mimics the absence of the pair production in this region due to the reduced magnetic field line curvature.
The default value is set to 30. This value should be at least one order of magnitude smaller than the maximum Lorentz factor acquired by the particles ($$\sim B_* \Omega^2$$, in normalized units).
`pp_gpair` is a parameter that defines the combined energy of the created particles in each pair production event. This energy is equally split between the secondary particles. The default value is set to 18. 
`pp_rmin` is a parameter that defines the minimum radius above which pairs can be created based on the selected threshold Lorentz factor. The default value is set to 1. 
`pp_rmax` is a parameter that defines the maximum radius below which pairs can be created based on the selected threshold Lorentz factor. The default value is set to 10. 

#### Available injection sites 
`inj_local` selects the spatial region where pairs are injected. The following options are available: 
- `"surface"` - Particles are injected at the stellar surface (lower radial boundary)
- `"volume"` - Particles are injected all over the simulation domain

Particles are injected through the selected injection mechanism (see `inj_type`) whenever the injection criterion is satisfied (see `inj_limit`). The default value is set to `"none"`. 
#### Injection associated parameters 
`inj_type` is a parameter that determines the charge density to inject depending on the preferred injection mechanism. The following injection mechanisms are available:
- `"fracnGJ"` - Particles are injected based on the local Goldreich-Julian (GJ) density - see e.g. ref. [1,4]
- `"fracnGJpole"` - Particles are injected based on the polar value of the GJ density 
- `"fracEB"` - Particles are injected based on the local value of the parallel electric field (parallel to the magnetic field) 
- `"fracEBteo"` - Particles are injected based on the local value of the theoretical vacuum parallel electric field (parallel to the magnetic field) 
- `"const"` - Particles are injected based on a specified constant density value 

The default value is set to `"none"`. 

`ks` determines the fraction of the injected charge density with respect to the selected reference (`inj_type`). For example, selecting `inj_type = "fracnGJ"` and `ks = 0.2`, we inject 20% of the local GJ density in each time step. The default value is set to 1. 
`inj_limit` is a parameter that defines the injection criterion, together with `k_limit` (see below). This parameter is used to determine a boolean that tells the code whether to inject new pairs or not in each cell. The following injection criteria are available: 
- `"nGJpole"` - Particles are injected if the total (electron+positron) local charge density is below `k_limit` times the GJ density on the pole
- `"altnGJpole"` - Particles are injected if the total (electron+positron) local number density is below `k_limit` times the GJ density on the pole
- `"nGJlocal"` - Particles are injected if the total (electron+positron) local charge density is below `k_limit` times the local GJ density
- `"altnGJlocal"` - Particles are injected if the total (electron+positron) local number density is below `k_limit` times the local GJ density
- `"fracEB"` - Particles are injected if the local value of the parallel electric field falls below `k_limit`
- `"fracEBteo"` - Particles are injected if the total (electron+positron) local number density falls below twice `k_limit` times the theoretical local vacuum parallel electric field value 
- `"sigma"` - Particles are injected if the local plasma magnetization falls below `sig_limit` (see below, only available for volume injection). 

The default value is set to `"none"`. 

`k_limit` is a parameter that represents the fraction of the physical quantity to compare in the injection criterion. The default value is set to 0. 

`sig_limit` is a parameter that defines the reference plasma magnetization. This parameter is only used in the injection criterion. The default value is set to 1000. 

`vs` is a parameter that defines the poloidal velocity with which particles are injected (only available for surface injection). This poloidal velocity is the velocity parallel to the poloidal magnetic field lines. When injecting in volume, particles are created with null momenta. The default value is set to 0. 

#### Available pusher types 
To push particles with general-relativistic effects, one needs to specify the push type. Here we present an example for the electrons: 
```
!-- species: electrons 
species 
{ 
  name = "electrons", 
  num_par_max = 100000, 
  rqm = -1.0, 
  num_par_x(1:2) = 1, 1, 
  add_tag = .false., 
  push_type = "rk4", 
} 
``` 

In the example above, we are selecting the 4th order Runge-Kutta method to push particles. The available types of pushers are: 
- `"standard"` - standard Boris pusher – for flat spacetime configurations (`rs = 0`) 
- `"rk4"` - Runge-Kutta 4th order pusher – for curved spacetime configurations (`rs > 1d-4`) 
- `"rk6"` - Runge-Kutta 6th order pusher – for curved spacetime configurations (`rs > 1d-4`) 

The default value is set to `"standard"`. 

## 5. References 
[1] A. Philippov, B. Cerutti, A. Tchekhovskoy, and Anatoly Spitkovsky, Astrophys. J. Lett. 815 L19 (2015)
[2] K. S. Thorne, R. H. Price, and D. A. Macdonald, Black Holes: The Membrane Paradigm (1986)
[3] L.D. Landau and E.M. Lifshitz, The Classical Theory of Fields, Pergamon (1975)
[4] B. Cerutti, A. Philippov, A., K. Parfrey, and A. Spitkovsky, Mon. Not. R. Astron. Soc. 448, 606 (2015) 
[5] F. Cruz, R. Torres, T. Grismayer, R. A. Fonseca, and L. O. Silva, in preparation
[6] L. Rezzolla, B. J. Ahmedov, and J. C. Miller, Mon. Not. R. Astron. Soc. 322, 723 (2001)