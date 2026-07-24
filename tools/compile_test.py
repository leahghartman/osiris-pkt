#!/usr/bin/env python

import subprocess
import sys, os, time
from optparse import OptionParser
import multiprocessing 
import itertools 

#-------------------------------------------------------------------------------
# Main script
#-------------------------------------------------------------------------------

# Parse command line options

parser = OptionParser( usage = "usage: %prog [options] [confs]" )
parser.add_option( "-v","--verbose", dest="verbose", action="store_true", default = False, \
				   help = "Verbose compilation" )
parser.add_option( "-d","--dims", dest="dims", metavar="DIMS" , default = '', \
				   help = "Comma separated values of dimensions to test (defaults to 1,2,3)" )
parser.add_option( "-p","--prec", dest="precs", metavar="PRECS" , default = '', \
				   help = "Comma separated values of precisions to test (options are SINGLE, DOUBLE. defaults to (SINGLE, DOUBLE)" )
parser.add_option( "-t","--type", dest="types", metavar="TYPE" , default = '', \
				   help = "Comma separated values of compilation type (defaults to production)" )
parser.add_option( "--docker_prefix", dest="docker_prefix", metavar="DOCKER_PREFIX" , default = '', \
				   help = "Script or command to put in front of all calls to make" )
parser.add_option( "--dry", dest = "dry", action="store_true", default = False,  \
				   help = "Script or command to put in front of all calls to make" )


(options, args) = parser.parse_args()

if ( len(args) == 0 ) :
	# No configuration(s) defined use default
	if ( sys.platform == 'darwin' ) :
		print('Using default configurations for darwin platform', file=sys.stderr)
		confs = ['macosx.gfortran', 'macosx.ifort64','macosx.pgi']
	elif ( sys.platform == 'linux2' ) :
		print('Using default configurations for linux platform', file=sys.stderr)
		confs = ['linux.gfortran', 'linux.ifort', 'linux.open64', 'linux.pgi'] 
	else :
		print('No configuration specified', file=sys.stderr)
		print('System:', sys.platform, file=sys.stderr)
		sys.exit('aborting...')

else :
	confs = args

# print( confs ) 

# Process dims
if ( options.dims == '' ) :
	dims = ['1','2','3']
else :
	dims = str(options.dims).split(',')

if ( options.precs == '' ) :
	precs = [ 'SINGLE', 'DOUBLE' ]
else :
	precs = str(options.dims).split(',')

# Process types
if ( options.types == '' ) :
	types = ['production', 'debug' ]
else :
	types = str(options.types).split(',')

docker_prefix = options.docker_prefix 

nproc = 8

# Print banner

print()
print("Testing compilation of osiris targets")
print()
print("Selected configurations : ")
print("-----------------------")
print("       Platform: ", confs )
print("     Dimensions: ", dims  )
print("          Types: ", types )
print("     Precisions: ", precs )
print()


report = []


def configure_and_compile( options ) : 

	conf, dim, _type, prec = options 

	# print('Testing compilation: ' + str( options ) )

	cwd = os.getcwd() + '/'

	# create a new directory and build there 
	builddir = './test/compilation-test/'
	builddir += '%s/%s/%s/%s/' % (conf, dim, _type, prec )
	os.makedirs( builddir, exist_ok = 1 )
	# os.chdir( builddir ) 

	# FNULL = open(os.devnull, 'w')

	cmd = '../' * 6  + 'configure -s %s -d %s -t %s' % (conf, dim, _type )
	# cmd = cwd + './configure -s %s -d %s -t %s' % (conf, dim, _type )
	# print( cmd ) 

	pipe = subprocess.Popen( cmd, shell = True, cwd = builddir,  
							stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
	
	out, err = pipe.communicate()

	if ( err ) :
		print('INFO: configure failed for ' + str( options ), file=sys.stderr )
		return 1

	# make
	t0 = time.time()
	
	cmd = 'make PRECISION=%s' % prec

	pipe = subprocess.Popen( cmd, shell = True, cwd = builddir, 
								# stdout= sys.stdout, 
								stdout=subprocess.PIPE, 
								stderr=subprocess.STDOUT  ) 
	
	out, err = pipe.communicate()

	# DEBUG 
	# out, err = 0, 0

	t1 = time.time()

	timediff = t1 - t0   # seconds 

	failure = err 

	if ( failure ) :
		print("[error] Compilation failed for configuration %s. \nError: %s" % ( str( options ), err ), file=sys.stderr)

	return ( failure, timediff ) 


compilation_test_failure = 0 
pool = multiprocessing.Pool( nproc )
results = []
failures = []
timediffs = [] 

# form iterable of cartesian product 
all_options = list( itertools.product( confs, dims, types, precs ) )

# # DEBUG : run one job 
# all_options = [ all_options[0] ]

# dispatch all jobs 
for options in all_options : 
	results.append( pool.apply_async( configure_and_compile, [options] ) )

pool.close() 

# collect and examine return codes 
for i in range( len( results ) ) : 
	failure, timediff = results[i].get() 
	failures.append( failure )
	timediffs.append( timediff )


if any( failures ) : 
	compilation_test_failure = 1 



print('') 
print('Compilation results:')
print('')

print('            Configuration    Dim         Type         Prec    Time [s]    Success')
print('-----------------------------------------------------------------')
for i in range( len( all_options ) ) :
	conf, dim, _type, prec = all_options[i]
	success = not failures[i]
	timediff = timediffs[i]
	print(" %24s %6s %12s %12s    %.1f         %d" % ( 
				conf, dim, _type, prec, timediff, success ) )


if compilation_test_failure : 
	sys.exit( 'ERROR: one or more compilations failed.' )

else : 
	print( '\n\nSUCCESS: all tests succeeded.')



