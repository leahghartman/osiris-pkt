import h5py, numpy as np
import sys



## PYTHON SCRIPT FOR CONVERTING IDL FORMATTED TRACKS TO HUMAN-READABLE FORMAT #### 

def main():
	# get args
	args = sys.argv[1:]

	# exit if no filename specified
	if(len(args) > 0):
		filename_in = args[0]
	else:
		print('No files specified')
		exit()


	try:
	    file_in = h5py.File(filename_in, 'r')
	except IOError:
	    print('cannot open ' + filename_in)
	    exit()

	# read data from file
	data = file_in['data'][:]
	itermap = file_in['itermap'][:]
	ntracks = file_in.attrs['NTRACKS'][0]
	niter = file_in.attrs['NITER'][0]
	quants = file_in.attrs['QUANTS'][:]
	file_in_attr_keys = file_in.attrs.keys()
	sim_attr_keys = file_in['SIMULATION'].attrs.keys()
	nquants = len(quants)

	# construct file out for new format
	filename_out = filename_in[:-3] + '-v2' + filename_in[-3:]
	file_out = h5py.File(filename_out,'w')

	# copy attrs from file_in
	for item in file_in_attr_keys:
		file_out.attrs[item] = file_in.attrs[item]
	for item in sim_attr_keys:
		file_out.attrs[item] = file_in['SIMULATION'].attrs[item]

	# first pass -- find total size of each track 
	#----------------------------------------#
	sizes = np.zeros(ntracks)

	itermapshape = itermap.shape
	for i in xrange(itermapshape[0]):
		part_number,npoints,nstart = itermap[i,:]
		sizes[part_number-1] += npoints

	#----------------------------------------#



	# initialize ordered data buffer
	#----------------------------------------#
	ordered_data = []
	for i in xrange(ntracks): 
		ordered_data.append(np.zeros((sizes[i],nquants)))
	#----------------------------------------#



	# assign tracks to ordered data from file_in data
	#----------------------------------------#
	track_indices = np.zeros(ntracks)
	data_index = 0

	for i in xrange(itermapshape[0]):
		part_number,npoints,nstart = itermap[i,:]
		track_index =  track_indices[part_number-1]


		ordered_data[part_number-1][track_index : track_index + npoints,0] \
		= nstart + np.arange(npoints) * niter

		ordered_data[part_number-1][track_index : track_index + npoints,1:] \
		= data[data_index:data_index + npoints ,:] 

		data_index += npoints
		track_indices[part_number-1] += npoints

	#----------------------------------------#

	# write to file out 
	for i in xrange(ntracks):
		group = file_out.create_group(str(i+1))
		for j in xrange(nquants):
			if(j==0):
				group.create_dataset(quants[j],data=np.array(ordered_data[i][:,j],dtype = int))
			else:
				group.create_dataset(quants[j],data=ordered_data[i][:,j])


	file_out.close()
	file_in.close()


if __name__ == "__main__":
	main()

