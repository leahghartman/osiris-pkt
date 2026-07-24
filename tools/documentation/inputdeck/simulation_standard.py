simulation_std = {
  
  'intro_text':  """ Hi There. Welcome to Osiris. This is the standard input description.
                  """
}


namelist_details = {
    '/osiris/read_sim_options/nl_simulation': { 
      'Label': 'Simulation' ,
      'detail-page-template': "standard",
      'summary-text':  """global simulation parameters""",
      'detail-page-header': """This section configures global simulation parameters and is optional. It accepts the following data:""",
      'detail-page-footer': """

      Here's an example of a simulation section defining a random number seed of 1237 and a reference density of 2.7e17 cm^-3.
      ```
        simulation 
        {
          random_seed = 1237, 
          n0 = 2.7e17,  ! [cm^-3]
        }
      ```""",
      'items': [

        ('algorithm', 
"""specifies the algorithm to be used in the simulation. 
   If set to "pgc" the ponderomotive guiding center solver (PGC) algorithm and 
   corresponding particle pusher are used. PGC is currently implemented in 
    2D geometries only. Please read node configuration and space configuration 
    sections for additional details on PGC usage. If set to "hd-hybrid" the hybrid algorithm is used."""),
        
        ('random_seed', 
"""specifies the seed used by the random number generator. 
   if not set it defaults to 0. If set the random number generator seed used by OSIRIS in each parallel 
   node will be random_seed + node_id (or just random_seed in serial runs)."""),

        ('omega_p0', 
"""specifies the reference simulation plasma frequency for the simulation in units of [rad/s]. 
    The user must set this value when doing simulations that require real (not simulation) units, 
    such as simulations involving ionization or collisions. Setting a reference density, n0, overrides 
    this value. See also the n0 option."""),

        ('n0', 
"""specifies the reference simulation plasma density for the simulation in units of [cm^-3]. 
    Setting this value overrides the omega_p0 option, calculating its value from the supplied value
    See also the omega_p0 option."""),

        ('gamma', 
"""specifies the relativistic factor of the simulation frame. This factor will transform all
   zpulses to the boosted frame, and will be used to automatically calculate the 'laserkdx' smooth coefficients."""),
        ('ndump_prof', 
""" specifies the frequency for profiling (timing) information dumps."""),
        ('wall_clock_limit', 
"""specifies a time limit to shut down the simulation, in a format "h:m:s". When 
   this wall clock time is reached the simulation will shut down. See also 'wall_clock_check'."""),
        ('wall_clock_check ', 
"""specifies the frequency at which to check if wall time limit has been reached. See also 'wall_clock_limit'."""),
        ('wall_clock_checkpoint', 
"""specifies if the code should dump checkpointing information when stopping due to wall clock limit.""")
      ]
    },

#------------------- Definitions for 'nl_node_conf'---------------------
#
 '/m_node_conf/read_nml_node_conf/nl_node_conf': { 
 'Label': 'Simulation' ,
 'detail-page-template': "standard",
 'summary-text':  """node configuration, and periodic boundary settings""",

 'detail-page-header': 
"""This section configures the parallel nodes and periodic boundary settings and must be present in the input file. 
It accepts the following data:""",

 'detail-page-footer': 
"""
Here's an example of a node_conf section for a 2D run, using 8 nodes in the x1 direction and 4 nodes in the x2 direction, 
and launching 2 threads per node. The run will use periodic boundaries in the x2 direction.
```
node_conf 
{
  node_number(1:2) = 8, 4,
  n_threads = 2,
  if_periodic(1:2) = .false., .true.,
}
```
""",
 'items': [
    ('node_number', 
"""
specifies the number of nodes to use in each direction for the simulation. The total number of nodes will be the product of the number 
of nodes for each direction. A single node run will be specified by setting all the items to 1. It is not required that the number of 
grid points on a given direction is evenly dividable by the number of nodes specified for that direction, but this will guaranty better 
load balancing. Note that only longitudinal partitions are allowed when algorithm = 'pgc' in the general simulation parameters section.
"""),
    ('n_threads', 
"""
specifies the number of threads per node to use. This option uses a shared memory algorithm inside each MPI node, that can significantly
improve load imbalance issues. To use this option the code must be compiled with OpenMP support. Setting this value to 1 will use the 
standard distributed memory algorithm only.
"""),
    ('if_periodic', 
"""
specifies if the boundary conditions for each direction will be periodic boundary conditions. If set, they will override any boundary 
conditions you specify for Electro-Magnetic fields or particle species.
"""),
    ('topology', 
"""
specifies how osiris should map simulation nodes to parallel nodes. Currently available options are:
- "mpi" - Use MPI functions (namely MPI_CART_CREATE) to generate the topology. MPI can use knowledge about network details to choose the best topology. This is the default.
- "old" - Use the old OSIRIS algorithm to generate the topology. This should only be used for testing purposes and comparing to old versions of the code.
- "bgq" - (available on BlueGene/Q systems only). This generates a topology based on the 5D torus topology of the BlueGene/Q system. 
          This should provide the best performance on these systems.
""")
 ]
 },





####################################################################
#                    Definitions for 'nl_grid'
####################################################################
 '/read_input_grid/nl_grid': { 
 'Label': 'Grid' ,
 'detail-page-template': "standard",
 'summary-text':  """grid and coordinate system""",

 'detail-page-header': 
"""This section configures the numerical grid and coordinate settings and must be present in the input file. 
It accepts the following data:""",

 'detail-page-footer': 
"""

""",
 'items': [
    ('nx_p ',
"""
"""),
    ('coordinates',
"""
"""),
    ('load_balance',
"""
"""),
    ('lb_type',
"""
"""),
    ('lb_gather',
"""
"""),
    ('n_dynamic',
"""
"""),
    ('start_load_balance',
"""
"""),
    ('max_imbalance',
"""
"""),
    ('cell_weight',
"""
"""),
    ('ndump_global_load',
"""
"""),
    ('ndump_node_load',
"""
"""),
    ('ndump_grid_load',
"""
"""),

 ]
 },

####################################################################
#           Defintions for nl_time_step
####################################################################
 '/m_time_step/read_nml_time_step/nl_time_step': {
 'Label':
"""
""",
 'detail-page-template':  'standard',
 'summary-text':
"""
""",
 'detail-page-header':
"""
""",
 'detail-page-footer':
"""
""",
 'items': [
    ('dt',
"""
"""),
    ('ndump',
"""
"""),
    ('dump_start',
"""
"""),
 ]
 },


####################################################################
#           Defintions for nl_restart
####################################################################
 '/m_restart/read_nml_restart/nl_restart': {
 'Label':
"""
""",
 'detail-page-template':  'standard',
 'summary-text':
"""
""",
 'detail-page-header':
"""
""",
 'detail-page-footer':
"""
""",
 'items': [
    ('ndump_fac',
"""
"""),
    ('if_restart',
"""
"""),
    ('if_remold',
"""
"""),
    ('ndump_time',
"""
"""),
    ('debug_iter',
"""
"""),
 ]
 },

####################################################################
#           Defintions for nl_space
####################################################################
 '/m_space/read_nml_space/nl_space': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('xmin',
""" 
"""),
    ('xmax',
""" 
"""),
    ('if_move',
""" 
"""),
    ('move_u',
""" 
"""),
 ]
 }, 

####################################################################
#           Defintions for nl_time
####################################################################
 '/m_time/read_nml_time/nl_time': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('tmin',
""" 
"""),
    ('tmax',
""" 
"""),
 ]
 },

####################################################################
#           Defintions for nl_el_mag_fld
####################################################################
 '/read_input_emf/nl_el_mag_fld': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('solver',
""" 
"""),
    ('k1',
""" 
"""),
    ('k2',
""" 
"""),
    ('smooth_type',
""" 
"""),
    ('smooth_niter',
""" 
"""),
    ('smooth_nmax',
""" 
"""),
    ('type_init_b',
""" 
"""),
    ('type_init_e',
""" 
"""),
    ('init_b0',
""" 
"""),
    ('init_e0',
""" 
"""),
    ('init_dipole_b_m',
""" 
"""),
    ('init_dipole_b_x0',
""" 
"""),
    ('init_dipole_b_r0',
""" 
"""),
    ('init_dipole_e_p',
""" 
"""),
    ('init_dipole_e_x0',
""" 
"""),
    ('init_dipole_e_r0',
""" 
"""),
    ('ext_fld',
""" 
"""),
    ('type_ext_b',
""" 
"""),
    ('type_ext_e',
""" 
"""),
    ('ext_b0',
""" 
"""),
    ('ext_e0',
""" 
"""),
    ('ext_b_mfunc',
""" 
"""),
    ('ext_e_mfunc',
""" 
"""),
    ('ext_dipole_b_m',
""" 
"""),
    ('ext_dipole_b_x0',
""" 
"""),
    ('ext_dipole_b_r0',
""" 
"""),
    ('ext_dipole_e_p',
""" 
"""),
    ('ext_dipole_e_x0',
""" 
"""),
    ('ext_dipole_e_r0',
""" 
"""),
    ('n_subcycle',
""" 
"""),
    ('marder_d',
""" 
"""),
    ('marder_n',
""" 
"""),
    ('fei_solver__n_coefficients',
""" 
"""),
    ('fei_solver__coefficients',
""" 
"""),
 ]
 },


####################################################################
#           Defintions for nl_emf_bound
####################################################################
 '/m_emf_bound/read_nml_emf_bound/nl_emf_bound': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('type',
""" 
"""),
    ('vpml_bnd_size',
""" 
"""),
    ('vpml_diffuse',
""" 
"""),
 ]
 },

####################################################################
#           Defintions for nl_smooth
####################################################################
 '/m_vdf_smooth/read_nml_smooth/nl_smooth': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('type',
""" 
"""),
    ('order',
""" 
"""),
    ('swfj',
""" 
"""),
    ('laserkdx_omega0',
""" 
"""),
    ('laserkdx_dir',
""" 
"""),
    ('digital_A',
""" 
"""),
    ('digital_fc',
""" 
"""),
 ]
 },

####################################################################
#           Defintions for nl_diag_emf
####################################################################
 '/m_emf_diag/read_input_diag_emf/nl_diag_emf': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('ndump_fac',
""" 
"""),
    ('ndump_fac_ave',
""" 
"""),
    ('ndump_fac_lineout',
""" 
"""),
    ('ndump_fac_ene_int',
""" 
"""),
    ('ndump_fac_charge_cons',
""" 
"""),
    ('prec',
""" 
"""),
    ('n_tavg',
""" 
"""),
    ('n_ave',
""" 
"""),
    ('reports',
""" 
"""),
 ]
 },

####################################################################
#           Defintions for nl_particles
####################################################################
 '/read_input_particles/nl_particles': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('num_species',
""" 
"""),
    ('num_cathode',
""" 
"""),
    ('num_neutral',
""" 
"""),
    ('num_neutral_mov_ions',
""" 
"""),
    ('low_jay_roundoff',
""" 
"""),
    ('ndump_fac',
""" 
"""),
    ('ndump_fac_ave',
""" 
"""),
    ('ndump_fac_lineout',
""" 
"""),
    ('n_ave',
""" 
"""),
    ('prec',
""" 
"""),
    ('reports',
""" 
"""),
    ('n_tavg',
""" 
"""),
    ('interpolation',
""" 
"""),
    ('grid_center',
""" 
"""),
    ('ndump_fac_ene',
""" 
"""),
 ]
 },

 ####################################################################
#           Defintions for nl_species
####################################################################
 '/read_input_species/nl_species': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('name',
""" 
"""),
    ('num_par_max',
""" 
"""),
    ('n_sort',
""" 
"""),
    ('rqm',
""" 
"""),
    ('q_real',
""" 
"""),
    ('num_par_x',
""" 
"""),
    ('tot_par_x',
""" 
"""),
    ('subcycle',
""" 
"""),
    ('push_type',
""" 
"""),
    ('push_start_time',
""" 
"""),
    ('num_pistons',
""" 
"""),
    ('add_tag',
""" 
"""),
    ('free_stream',
""" 
"""),
    ('init_fields',
""" 
"""),
    ('if_collide',
""" 
"""),
    ('if_like_collide',
""" 
"""),
    ('init_type',
""" 
"""),
 ]
 },

####################################################################
#           Defintions for nl_udist
####################################################################
 '/m_species_udist/read_nml_spec_udist/nl_udist': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('uth_type',
""" 
"""),
    ('use_spatial_uth',
""" 
"""),
    ('use_spatial_ufl',
""" 
"""),
    ('uth',
""" 
"""),
    ('ufl',
""" 
"""),
    ('spatial_uth',
""" 
"""),
    ('spatial_ufl',
""" 
"""),
    ('relmax_T',
""" 
"""),
    ('relmax_umax',
""" 
"""),
    ('use_classical_uadd',
""" 
"""),
    ('n_accelerate',
""" 
"""),
    ('n_q_incr',
""" 
"""),
    ('n_accelerate_type',
""" 
"""),
    ('n_q_incr_type',
""" 
"""),
 ]
 },

####################################################################
#           Defintions for nl_spe_bound
####################################################################
 '/m_species_boundary/read_nml_spe_bound/nl_spe_bound': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('type',
""" 
"""),
    ('uth_bnd',
""" 
"""),
    ('ufl_bnd',
""" 
"""),
    ('thermal_type',
""" 
"""),
 ]
 },

####################################################################
#           Defintions for nl_piston
####################################################################
 '/m_piston/read_nml_pst/nl_piston': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('dim',
""" 
"""),
    ('updown',
""" 
"""),
    ('v',
""" 
"""),
    ('start_pos',
""" 
"""),
    ('start_time',
""" 
"""),
    ('stop_time',
""" 
"""),
    ('opacity_factor',
""" 
"""),
    ('profile_type',
""" 
"""),
    ('profile_expression',
""" 
"""),
 ]
 },

####################################################################
#           Defintions for nl_diag_species
####################################################################
 '/read_input_diag_species/nl_diag_species': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('ndump_fac_pha',
""" 
"""),
    ('ndump_fac_pha_tavg',
""" 
"""),
    ('ndump_fac_lineout',
""" 
"""),
    ('ndump_fac_ene',
""" 
"""),
    ('ndump_fac_heatflux',
""" 
"""),
    ('ndump_fac_temp',
""" 
"""),
    ('ndump_fac_raw',
""" 
"""),
    ('ndump_fac',
""" 
"""),
    ('ndump_fac_ave',
""" 
"""),
    ('n_ave',
""" 
"""),
    ('prec',
""" 
"""),
    ('ps_xmin',
""" 
"""),
    ('ps_xmax',
""" 
"""),
    ('ps_pmin',
""" 
"""),
    ('ps_pmax',
""" 
"""),
    ('ps_lmin',
""" 
"""),
    ('ps_lmax',
""" 
"""),
    ('if_ps_p_auto',
""" 
"""),
    ('if_ps_l_auto',
""" 
"""),
    ('ps_gammamin',
""" 
"""),
    ('ps_gammamax',
""" 
"""),
    ('if_ps_gamma_log',
""" 
"""),
    ('if_ps_gamma_auto',
""" 
"""),
    ('ps_nx',
""" 
"""),
    ('ps_nx_3D',
""" 
"""),
    ('ps_np',
""" 
"""),
    ('ps_np_3D',
""" 
"""),
    ('ps_nl',
""" 
"""),
    ('ps_nl_3D',
""" 
"""),
    ('ps_ngamma',
""" 
"""),
    ('raw_gamma_limit',
""" 
"""),
    ('raw_fraction',
""" 
"""),
    ('raw_math_expr',
""" 
"""),
    ('n_ene_bins',
""" 
"""),
    ('ene_bins',
""" 
"""),
    ('ndump_fac_tracks',
""" 
"""),
    ('n_start_tracks',
""" 
"""),
    ('niter_tracks',
""" 
"""),
    ('file_tags',
""" 
"""),
    ('ifdmp_tracks_efl',
""" 
"""),
    ('ifdmp_tracks_bfl',
""" 
"""),
    ('ifdmp_tracks_psi',
""" 
"""),
    ('phasespaces',
""" 
"""),
    ('pha_ene_bin',
""" 
"""),
    ('pha_cell_avg',
""" 
"""),
    ('pha_time_avg',
""" 
"""),
    ('n_tavg',
""" 
"""),
    ('reports',
""" 
"""),
    ('rep_cell_avg',
""" 
"""),
    ('rep_udist',
""" 
"""),
 ]
 },

####################################################################
#           Defintions for nl_cathode
####################################################################
 '/m_cathode/read_input_cathode/nl_cathode': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('dir',
""" 
"""),
    ('wall',
""" 
"""),
    ('t_start',
""" 
"""),
    ('t_rise',
""" 
"""),
    ('t_fall',
""" 
"""),
    ('t_flat',
""" 
"""),
    ('mov_inj',
""" 
"""),
    ('prof_type',
""" 
"""),
    ('density',
""" 
"""),
    ('center',
""" 
"""),
    ('gauss_width',
""" 
"""),
    ('gauss_w0',
""" 
"""),
    ('tr_u_expression',
""" 
"""),
    ('channel_bottom',
""" 
"""),
    ('channel_r0',
""" 
"""),
    ('channel_depth',
""" 
"""),
    ('channel_size',
""" 
"""),
    ('channel_wall',
""" 
"""),
    ('deposit_current',
""" 
"""),
    ('den_min',
""" 
"""),
 ]
 },

####################################################################
#           Defintions for nl_neutral
####################################################################
 '/m_neutral/read_input_neutral/nl_neutral': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('name',
""" 
"""),
    ('neutral_gas',
""" 
"""),
    ('ion_param',
""" 
"""),
    ('den_min',
""" 
"""),
    ('e_min',
""" 
"""),
    ('multi_max',
""" 
"""),
    ('multi_min',
""" 
"""),
    ('if_tunnel',
""" 
"""),
    ('if_impact',
""" 
"""),
    ('inject_line',
""" 
"""),
 ]
 },

####################################################################
#           Defintions for nl_neutral_mov_ions
####################################################################
 '/m_neutral/read_input_neutral/nl_neutral_mov_ions': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('name',
""" 
"""),
    ('neutral_gas',
""" 
"""),
    ('ion_param',
""" 
"""),
    ('den_min',
""" 
"""),
    ('e_min',
""" 
"""),
    ('multi_max',
""" 
"""),
    ('multi_min',
""" 
"""),
    ('if_tunnel',
""" 
"""),
    ('if_impact',
""" 
"""),
    ('inject_line',
""" 
"""),
 ]
 },

####################################################################
#           Defintions for nl_diag_neutral
####################################################################
 '/m_diag_neutral/read_nml_diag_neutral/nl_diag_neutral': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('ndump_fac',
""" 
"""),
    ('ndump_fac_ave',
""" 
"""),
    ('ndump_fac_lineout',
""" 
"""),
    ('prec',
""" 
"""),
    ('n_ave',
""" 
"""),
    ('n_tavg',
""" 
"""),
    ('reports',
""" 
"""),
 ]
 },

####################################################################
#           Defintions for nl_collisions
####################################################################
 '/m_species_collisions/read_nml_coll/nl_collisions': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('n_collide',
""" 
"""),
    ('nx_collision_cells',
""" 
"""),
    ('coulomb_logarithm_automatic',
""" 
"""),
    ('coulomb_logarithm_value',
""" 
"""),
    ('collision_model',
""" 
"""),
    ('cross_section_correction',
""" 
"""),
    ('timestep_cross_corr',
""" 
"""),
    ('frame_correction',
""" 
"""),
    ('root_finder',
""" 
"""),
    ('rng_seed',
""" 
"""),
    ('perez_low_temp_correction',
""" 
"""),
 ]
 },


####################################################################
#           Defintions for nl_zpulse
####################################################################
 '/m_zpulse_std/read_input_zpulse/nl_zpulse': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('if_launch',
""" 
"""),
    ('b_type',
""" 
"""),
    ('a0',
""" 
"""),
    ('omega0',
""" 
"""),
    ('phase',
""" 
"""),
    ('pol_type',
""" 
"""),
    ('pol',
""" 
"""),
    ('propagation',
""" 
"""),
    ('direction',
""" 
"""),
    ('chirp_order',
""" 
"""),
    ('chirp_coefs',
""" 
"""),
    ('lon_fwhm',
""" 
"""),
    ('lon_type',
""" 
"""),
    ('lon_start',
""" 
"""),
    ('lon_rise',
""" 
"""),
    ('lon_flat',
""" 
"""),
    ('lon_fall',
""" 
"""),
    ('lon_duration',
""" 
"""),
    ('lon_x0',
""" 
"""),
    ('lon_range',
""" 
"""),
    ('lon_math_func',
""" 
"""),
    ('lon_tilt',
""" 
"""),
    ('per_type',
""" 
"""),
    ('per_center',
""" 
"""),
    ('per_w0',
""" 
"""),
    ('per_fwhm',
""" 
"""),
    ('per_focus',
""" 
"""),
    ('per_chirp_order',
""" 
"""),
    ('per_chirp_coefs',
""" 
"""),
    ('per_w0_asym',
""" 
"""),
    ('per_fwhm_asym',
""" 
"""),
    ('per_focus_asym',
""" 
"""),
    ('per_asym_trans',
""" 
"""),
    ('per_n',
""" 
"""),
    ('per_0clip',
""" 
"""),
    ('per_kt',
""" 
"""),
    ('launch_time',
""" 
"""),
    ('no_div_corr',
""" 
"""),
    ('per_tem_mode',
""" 
"""),
 ]
 },


####################################################################
#           Defintions for nl_diag_current
####################################################################
 '/m_current_diag/read_input_diag_current/nl_diag_current': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('ndump_fac',
""" 
"""),
    ('ndump_fac_ave',
""" 
"""),
    ('ndump_fac_lineout',
""" 
"""),
    ('prec',
""" 
"""),
    ('n_ave',
""" 
"""),
    ('n_tavg',
""" 
"""),
    ('reports',
""" 
"""),
 ]
 },


####################################################################
#           Defintions for nl_antenna_array
####################################################################
 '/m_antenna_array/read_nml_antenna_array/nl_antenna_array': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('n_antenna',
""" 
"""),
 ]
 },
 

####################################################################
#           Defintions for nl_antenna
####################################################################
 '/m_antenna/read_nml_antenna/nl_antenna': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('a0',
""" 
"""),
    ('t_rise',
""" 
"""),
    ('t_flat',
""" 
"""),
    ('t_fall',
""" 
"""),
    ('t_offset',
""" 
"""),
    ('omega0',
""" 
"""),
    ('focus',
""" 
"""),
    ('x0',
""" 
"""),
    ('y0',
""" 
"""),
    ('rad_x',
""" 
"""),
    ('rad_y',
""" 
"""),
    ('ant_type',
""" 
"""),
    ('spin',
""" 
"""),
    ('pol',
""" 
"""),
    ('delay',
""" 
"""),
    ('side',
""" 
"""),
    ('tilt',
""" 
"""),
    ('phase',
""" 
"""),
    ('direction',
""" 
"""),
    ('power',
""" 
"""),
    ('eps',
""" 
"""),
    ('chirp',
""" 
"""),
    ('jitterw',
""" 
"""),
    ('jitteromg',
""" 
"""),
 ]
 },

 '/m_psource_btfel/read_input_betatronfel/nl_betatron_fel': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('gamma_0',
""" 
"""),
    ('K_number',
""" 
"""),
    ('DE_E',
""" 
"""),
    ('DK_K',
""" 
"""),
    ('rotation_center',
""" 
"""),
    ('init_random',
""" 
"""),
    ('in_boosted_frame',
""" 
"""),
 ]
 },


 ####################################################################
#           Defintions for nl_num_ene
####################################################################
 '/m_cross/read_nml_cross/nl_num_ene': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('num_ene',
""" 
"""),
 ]
 },
####################################################################
#           Defintions for nl_cross
####################################################################
 '/m_cross/read_nml_cross/nl_cross': {
 'Label': 
""" 
""",
 'detail-page-template':  'standard', 
 'summary-text':
""" 
""",
 'detail-page-header':
""" 
""",
 'detail-page-footer':
""" 
""",
 'items': [
    ('cross',
""" 
"""),
    ('ene',
""" 
"""),
 ]
 },
 
####################################################################              
#           Defintions for nl_profile                                             
####################################################################              
 '/m_psource_std/read_input_std/nl_profile': {                                    
 'Label':                                                                         
"""                                                                               
""",                                                                              
 'detail-page-template':  'standard',                                             
 'summary-text':                                                                  
"""                                                                               
""",                                                                              
 'detail-page-header':                                                            
"""                                                                               
""",                                                                              
 'detail-page-footer':                                                            
"""                                                                               
""",                                                                              
 'items': [                                                                       
    ('profile_type',                                                              
"""                                                                               
"""),                                                                             
    ('den_min',                                                                   
"""                                                                               
"""),                                                                             
    ('density',                                                                   
"""                                                                               
"""),                                                                             
    ('num_x',                                                                     
"""                                                                               
"""),                                                                             
    ('x',                                                                         
"""                                                                               
"""),                                                                             
    ('fx',                                                                        
"""                                                                               
"""),                                                                             
    ('gauss_n_sigma',                                                             
"""                                                                               
"""),                                                                             
    ('gauss_range',                                                               
"""                                                                               
"""),                                                                             
    ('gauss_center',                                                              
"""                                                                               
"""),                                                                             
    ('gauss_sigma',                                                               
"""                                                                               
"""),                                                                             
    ('channel_dir',                                                               
"""                                                                               
"""),                                                                             
    ('channel_r0',                                                                
"""                                                                               
"""),                                                                             
    ('channel_depth',                                                             
"""                                                                               
"""),                                                                             
    ('channel_size',                                                              
"""                                                                               
"""),                                                                             
    ('channel_center',                                                            
"""                                                                               
"""),                                                                             
    ('channel_wall',                                                              
"""                                                                               
"""),                                                                             
    ('channel_pos',                                                               
"""                                                                               
"""),                                                                             
    ('channel_bottom',                                                            
"""                                                                               
"""),                                                                             
    ('sphere_center',                                                             
"""                                                                               
"""),                                                                             
    ('sphere_radius',                                                             
"""                                                                               
"""),                                                                             
    ('math_func_expr',                                                            
"""                                                                               
"""),                                                                             
    ('sample_rate',                                                               
"""                                                                               
"""),                                                                             
 ]                                                                                
 },                                                                               
# --- close 'namelist_details' dictionary                                         
}



