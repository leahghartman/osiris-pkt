import h5py
import numpy as np
import scipy.constants as const

def osirisLoadDensity(  dump , n0 = 1. , species = 'plasma', path = './' ):
	import h5py
	import numpy as np
	

	n = h5py.File( path + 'MS/DENSITY/%s/charge/charge-%s-%.6d.h5' % (species, species, dump) , 'r')['charge']

	return np.asarray( n )
		

def osirisLoadEField(  dump , component = '1' , e0 = 1. , modes = [0,1,2] , path = './' ):
	import h5py
	import numpy as np
	

	E = h5py.File( path + 'MS/FLD/e%s/e%s-%.6d.h5' % (component,component,dump) , 'r' )['e%s' % component]

			
	return np.asarray( E )
	

def osirisLoadBField(  dump , component = '1' , b0 = 1. , modes = [0,1,2,3,4] , path = './' ):
	import h5py
	import numpy as np
	

	B = h5py.File( path + 'MS/FLD/b%s/b%s-%.6d.h5' % (component,component,dump) , 'r' )['b%s' % component]

			
	return np.asarray( B )
	

def osirisLoadPhaseSpace(  dump , species = 'plasma', phasespace = 'p1x1', path = './' ):
	import h5py
	import numpy as np
	    
	out = h5py.File( path + 'MS/PHA/%s/%s/%s-%s-%.6d.h5' % (phasespace,species,phasespace,species,dump) , 'r')
	f = np.array(out[phasespace])
	    
	axes = []
	N_axes = len(out['AXIS'].keys())
	for dim in range(N_axes):
		nax = np.shape(f)[N_axes - (dim+1)]
		axmin, axmax = out['AXIS']['AXIS%d'%(dim+1)]
		axes.append(np.linspace(axmin, axmax, nax)) 
	     
	return axes, f
		
	
    
def osirisLoadGrid(  dump  , x0 = 1.  , species = 'plasma', path = './' ):
	import h5py
	import numpy as np

	f = h5py.File( path + 'MS/DENSITY/%s/charge/charge-%s-%.6d.h5' % (species, species,dump) )
	
	
	Nz, = f['charge'].shape
	zmin , zmax = f['AXIS']['AXIS1'][:]
	
	
	dz = ( zmax - zmin ) / Nz
	
	z = x0 * np.linspace( zmin  , zmax  , Nz )
	
	
	return z
	
def osirisGridSpacing(  dump  , x0 = 1.  , species = 'plasma', path = './' ):
	import h5py
	import numpy as np

	f = h5py.File( path + 'MS/DENSITY/%s/charge/charge-%s-%.6d.h5' % (species, species,dump) )
	
	
	Nz,   = f['charge'].shape
	
	zmin , zmax = f['AXIS']['AXIS1'][:]
	
	dz = ( zmax - zmin ) / Nz
	
	return dz





#########################################################################


def FieldSlice( field , r , theta ):
	import numpy as np


	f1 = fbpicFieldReconstructionFromNodes( field , theta         )
	f2 = fbpicFieldReconstructionFromNodes( field , theta + np.pi )



	x  = np.concatenate( [ - r[::-1]  , r  ] , axis = 0 )
	f  = np.concatenate( [  f2[::-1]  , f1 ] , axis = 0 )

	return x , f


def E_xy_Slices( E_r , E_t , r , theta ):
	import numpy as np

	E_r_pos = fbpicFieldReconstructionFromNodes( E_r , theta         )
	E_r_neg = fbpicFieldReconstructionFromNodes( E_r , theta + np.pi )

	E_t_pos = fbpicFieldReconstructionFromNodes( E_t , theta         )
	E_t_neg = fbpicFieldReconstructionFromNodes( E_t , theta + np.pi )


	E_x_pos = E_r_pos * np.cos(theta) - E_t_pos * np.sin(theta)
	E_y_pos = E_r_pos * np.sin(theta) + E_t_pos * np.cos(theta)

	E_x_neg = E_r_neg * np.cos( theta + np.pi ) - E_t_neg * np.sin( theta + np.pi )
	E_y_neg = E_r_neg * np.sin( theta + np.pi ) + E_t_neg * np.cos( theta + np.pi )

	x    = np.concatenate( [ - r[::-1]       , r       ] , axis = 0 )
	E_x  = np.concatenate( [  E_x_neg[::-1]  , E_x_pos ] , axis = 0 )
	E_y  = np.concatenate( [  E_y_neg[::-1]  , E_y_pos ] , axis = 0 )

	return x , E_x , E_y

def B_xy_Slices( B_r , B_t , r , theta ):
	import numpy as np

	B_r_pos = fbpicFieldReconstructionFromNodes( B_r , theta         )
	B_r_neg = fbpicFieldReconstructionFromNodes( B_r , theta + np.pi )

	B_t_pos = fbpicFieldReconstructionFromNodes( B_t , theta         )
	B_t_neg = fbpicFieldReconstructionFromNodes( B_t , theta + np.pi )


	B_x_pos = B_r_pos * np.cos(theta) - B_t_pos * np.sin(theta)
	B_y_pos = B_r_pos * np.sin(theta) + B_t_pos * np.cos(theta)

	B_x_neg = B_r_neg * np.cos( theta + np.pi ) - B_t_neg * np.sin( theta + np.pi )
	B_y_neg = B_r_neg * np.sin( theta + np.pi ) + B_t_neg * np.cos( theta + np.pi )

	x    = np.concatenate( [ - r[::-1]       , r       ] , axis = 0 )
	B_x  = np.concatenate( [  B_x_neg[::-1]  , B_x_pos ] , axis = 0 )
	B_y  = np.concatenate( [  B_y_neg[::-1]  , B_y_pos ] , axis = 0 )

	return x , B_x , B_y






##########################################################################


def fbpicToVTKdirect( theta , r , z , diags = {} , filename = 'para'  ):
	import numpy as np
	from evtk.hl import gridToVTK

	"""
	Save fbpic grid diagnostics

	INPUTS:
	=======
	theta = linear array of azimuthal grid points , len(theta) = Ntheta
	r     = linear array of radial grid points    , len(r)     = Nr
	z     = linear array of axial grid points     , len(z)     = Nz

	filename = string 

	diags: dictionary of diagnostics
		keys   -- will appear as the key in paraview
		values -- <HDF5 dataset "r": shape ( 2*Nmodes + 1 , Nr , Nz )
	
		Nz and Nr have to be consistent with the sizes of r and z


	example:

	diags = { 'rho (C/cm^3)' : - density[...,::10] / 1e6 , # C/cm^3
			  'Ez (GV/m)' :   Ez[...,::10] /1e9            # GV/m
											  }

	GridReconstruction will be applied to each diagnostic, then the
	diagnostics are saved as a VTK structured array

	
	With this design it is difficult to save Ex and Ey because they need
	to be constructed from Er and Et after the GridReconstruction

	"""    

	thetagrid,rgrid,zgrid =  np.meshgrid(   theta , r ,  z , indexing =  'ij' )

	xgrid = rgrid * np.cos( thetagrid )
	ygrid = rgrid * np.sin( thetagrid )


	# print 'xgrid.shape:' , xgrid.shape
	# print 'ygrid.shape:' , ygrid.shape
	# print 'zgrid.shape:' , zgrid.shape
		

	reconstructed_diags = {}


	for diagname,diagdata in diags.items():

		diagdata3d = fbpicFieldReconstructionFromNodes( field = diagdata , theta = thetagrid )
		# print '%s.shape' % diagname , diagdata3d.shape

		reconstructed_diags[ diagname ] = diagdata3d


	gridToVTK( filename , xgrid, 
						  ygrid, 
						  zgrid, 
						  pointData = reconstructed_diags )


	return True




def osiris_density_slice(dump, theta = 0, n0=1., species = 'plasma', modes = [0,1,2,3,4], path = './' ):
	import numpy as np
	number_density_modes  = osirisLoadDensity( dump = dump , n0 = n0 , species = species, modes = modes, path = path ) # in 1/cm^3
	number_density_osiris1 = osirisFieldReconstructionFromNodes( number_density_modes , theta = theta ) # in 1/cm^3
	number_density_osiris2 = osirisFieldReconstructionFromNodes( number_density_modes , theta = (theta + np.pi)) # in 1/cm^3
	return  np.vstack((number_density_osiris1[::-1],number_density_osiris2))



def osirisLoadGrid2r(  dump  , x0 = 1.  , species = 'plasma',  path = './' ):
	import h5py
	import numpy as np

	f = h5py.File( path + 'MS/DENSITY/%s/MODE-0-RE/charge_cyl_m/charge_cyl_m-%s-0-re-%.6d.h5' % (species, species, dump) )
	
	
	Nr   , Nz   = f['charge_cyl_m'].shape
	
	zmin , zmax = f['AXIS']['AXIS1'][:]
	rmin , rmax = f['AXIS']['AXIS2'][:]
	
	dz = ( zmax - zmin ) / Nz
	dr = ( rmax - rmin ) / Nr
		
	z = x0 * np.linspace( zmin  , zmax  , Nz ) #+ 0.5 * dz
	r = x0 * np.linspace( -rmax  , rmax  , 2 * Nr ) 
	
	return r , z

def osiris_Ez_slice(dump, theta = 0, e0 = 1.0, modes = [0,1,2],  species = 'plasma',path = './', x0 = 1, ):
	Ez_modes     = osirisLoadEField( dump, component = '1' , e0 = e0 , modes = modes, path = path )  
	r , _        = osirisLoadGrid(   dump = dump , x0 = x0,  species = species, path = path)
	_ , Ez_slice = FieldSlice( Ez_modes , r , theta = theta )
	return Ez_slice

def osiris_Bz_slice(dump, theta = 0, b0 = 1.0, modes = [0,1,2],  species = 'plasma',path = './', x0 = 1, ):
	Bz_modes     = osirisLoadBField( dump, component = '1' , b0 = b0 , modes = modes, path = path )  
	r , _        = osirisLoadGrid(   dump = dump , x0 = x0,  species = species, path = path)
	_ , Bz_slice = FieldSlice( Bz_modes , r , theta = theta )
	return Bz_slice

def osiris_Exy_slice(dump, theta = 0, e0 = 1.0, modes = [0,1,2],  species = 'plasma', path = './', x0 = 1, ):
	Er_modes = osirisLoadEField( dump,  component = '2' , e0 = e0 , modes = modes,  path = path )  
	Et_modes = osirisLoadEField( dump , component = '3' , e0 = e0 , modes = modes , path = path )
	r, _     = osirisLoadGrid( dump = dump , x0 = x0,  species = species, path = path)
	_, Ex, Ey= E_xy_Slices( Er_modes , Et_modes , r , theta )
	return Ex, Ey

def osiris_Bxy_slice(dump, theta = 0, b0 = 1.0, modes = [0,1,2,3,4],  species = 'plasma', path = './', x0 = 1, ):
	Br_modes = osirisLoadBField( dump,  component = '2' , b0 = b0 , modes = modes,  path = path )  
	Bt_modes = osirisLoadBField( dump , component = '3' , b0 = b0 , modes = modes , path = path )
	r, _     = osirisLoadGrid( dump = dump , x0 = x0,  species = species, path = path)
	_, Bx, By= B_xy_Slices( Br_modes , Bt_modes , r , theta )
	return Bx, By

def osiris_particles(dump, particles = 'plasma',path = './', x0 = 1 ):
	import h5py
	import numpy as np
	f = h5py.File( path + 'MS_pkt/RAW/%s/RAW-%s-%.6d.h5' % (particles,particles,dump) )
	ene = f['ene'][()]
	p1 = f['p1'][()]
	p2 = f['p2'][()]
	p3 = f['p3'][()]
	wgt= f['q'][()]
	z  = f['x1'][()] * x0
	return np.asarray([z,p1,p2,p3,ene,wgt])




def fb_diag(path, dump, deck, field = 'E', density = 'rho_e1' ):
	import h5py
	data = h5py.File(path+ 'data%08d.h5' % dump, 'r' )
	density     = data['data']['%d' % dump ]['fields']['%s'%density]
	fz          = data['data']['%d' % dump ]['fields']['%s'%field]['z']

	r = np.linspace(0, deck.rmax,deck.Nr)*1e6
	z = (np.linspace(deck.zmin, deck.zmax, deck.Nz)+ const.c * deck.dt * dump) / 1e-6
	
	x, density_rec = FieldSlice(density , r , 0) 
	density_rec = -density_rec/const.elementary_charge / 1e24
	_, fz_rec = FieldSlice(fz , r , theta = 0) 
	fz_rec = fz_rec/1e9
	
	return x, z, density_rec, fz_rec

def fb_diag_xy(path, dump, deck, field = 'E', density = 'rho_e1' ):
	import h5py
	data = h5py.File(path+ 'data%08d.h5' % dump, 'r' )
	density     = data['data']['%d' % dump ]['fields']['%s'%density]
	fr          = data['data']['%d' % dump ]['fields']['%s'%field]['r']
	ft          = data['data']['%d' % dump ]['fields']['%s'%field]['t']
	
	r = np.linspace(0, deck.rmax,deck.Nr)*1e6
	
	x, fx, fy = E_xy_Slices( fr , ft , r , theta = 0 )
	
	return x, fx / 4e12, fy / 4e12



def fb_diag_exy(path, dump, deck, density = 'e1' ):
    data = h5py.File(path+ 'data%08d.h5' % dump, 'r' )
    density     = data['data']['%d' % dump ]['fields']['%s'%density]
    Er          = data['data']['%d' % dump ]['fields']['E']['r']
    Et          = data['data']['%d' % dump ]['fields']['E']['t']
    
    r = np.linspace(0, deck.rmax,deck.Nr)*1e6
    
    x, Ex, Ey = E_xy_Slices( Er , Et , r , theta = 0 )
    
    return x, Ex / 4e12, Ey / 4e12




