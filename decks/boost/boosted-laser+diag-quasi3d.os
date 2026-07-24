simulation
{
	gamma = 10.0,
	algorithm = "quasi-3D",
}

node_conf
{
	node_number(1:2) = 2, 2,
	n_threads = 1,
	if_periodic(1:2) = .false., .false.,
}

grid
{
	nx_p(1:2) = 312, 30,
	coordinates = "cylindrical",
	n_cyl_modes = 1,
}

time_step
{
	dt = 0.09974937185533099,
	ndump = 50,
}

restart
{
	ndump_fac = 15,
	if_restart = .false.,
	if_remold = .true.,
}

space
{
	xmin(1:2) = -100.0, 0.0,
	xmax(1:2) = 24.2, 12.0,
	if_move = .false., .false.,
}

time
{
	tmin = 0.0,
	tmax = 112,
}

el_mag_fld
{
	solver = "custom",
}

emf_bound
{
	type(1:2,1) = "lindman", "lindman",
	type(1:2,2) = "axial", "open",
}

emf_solver
{
	type = "dual",
	solver_ord = 2,
	n_coef = 16,
	weight_n = 10,
	weight_w = 0.3,
	filter_limit = 0.6,
	filter_width = 0.1,
	n_damp_cell = 10,
	filter_current = .true.,
	correct_current = .true.,
}

diag_emf
{
	ndump_fac = 1,
	ndump_fac_ave = 0,
	ndump_fac_lineout = 0,
	reports = "e1_cyl_m", "e2_cyl_m", "e3_cyl_m", "s1_int_cyl_m",
}

boosted_diag
{
	ndump_fac = 1,
	boost_tmin = 0.0,
	boost_tmax = 240.0,
	boost_dt = 6.0,
	boost_xmin = -10.0,
	boost_xmax = 2.0,
	boost_nx = 600,
	boost_if_move = .true.,
	boost_move_vel = 1.0,
	boosted_reports = "e1_boost_cyl_m", "e2_boost_cyl_m", "e3_boost_cyl_m", "b1_boost_cyl_m", "b2_boost_cyl_m", "b3_boost_cyl_m",
}

particles
{
	num_species = 0,
	interpolation = "linear",
}

zpulse_mov_wall
{
	a0 = 4.0,
	if_launch = .true.,
	omega0 = 10,
	pol_type = 0,
	pol = 0.0,
	wall_pos = -6.782957286162508,
	wall_vel = -0.99498743710662,
	ncells = 17,
	tenv_type = "polynomial",
	tenv_rise = 4.0,
	tenv_fall = 4.0,
	per_type = "gaussian",
	per_w0 = 4.0,
	xi0 = 0,
	per_focus = 0,
}


current
{
}

smooth
{
}

