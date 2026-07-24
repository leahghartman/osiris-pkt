# CUDA algorithm

## Overview

The CUDA module enables one to run OSIRIS on NVIDIA GPU-accelerated architectures. It implements a GPU parallelized particle pusher and current deposit in CUDA, while relying on the existing CPU routines for the field solve, particle communication, particle boundary conditions, field boundary conditions, and particle initialization. Particle sorting is also done on the GPU so that particles can live in device memory most of the time and expensive copy operations between device and host can be kept to a minimum.

CUDA is implemented as a simulation mode which inherits from the tiles simulation mode. For this reason, some of the features and parameters used here will be new to a user of "standard" OSIRIS. See `source/tiles/README.md` for a full description of the tiles simulation mode. 

In the remainder of the "Overview" section, some of the inner-workings of the CUDA algorithm, along with their associated input deck parameters, are described. For most cases, it is not necessary to understand this stuff, and one can simply use the default values of these parameters. However, certain problems may cause the code to produce an error message and exit the simulation. In these cases, some knowledge of the inner-workings of the code can be useful for diagnosing and fixing the problem.

To efficiently store particles in GPU memory, each MPI rank stores particles on its respective GPU in "chunks" of memory which store data for a specified number of particles. Each rank allocates most of the device memory at initialization as chunks, which can then be freely grabbed from and returned to the global chunk pool as necessary with little overhead. Each tile maintains a "pad" of enough extra chunks of memory such that buffer overflow during communication between tiles is unlikely, but few enough such that each tile nearly saturates its memory capacity (i.e. wastes little memory). Ideally, this results in maximal device memory utilization: each tile has just the right amount of memory, and there is nearly zero overhead to effectively resizing its memory by grabbing or returning chunks.

In practice, the chunk pool is controlled by two parameters: `chunk_size`, which is the number of particles each chunk can store, and `chunk_pool_size`, which is the total number of chunks allocated by each MPI rank. Using a small chunk size results in greater flexibility for memory management, but additional overhead for grabbing and returning chunks from the pool. We have found that the performance penalty from the latter is neglibigle; thus, we recommend using the smallest possible chunk size of 512. On the other hand, `chunk_pool_size` must be chosen such that `chunk_size` * `chunk_pool_size` $\geq$ N, where N is the number of particles on each GPU. Note that with dynamic load balancing, N can be assumed to be approximately constant. In practice, the default values should be appropriate for most applications. This sets the `chunk_size = 512` and the `chunk_pool_size` such that chunks take up most of the available memory on the GPU (75%).

As mentioned above, in addition to holding enough chunks to fit all its particles, each tile holds some additional chunks for buffering particles during the particle sort. Specifically, each tile holds `nchunks_pad` chunks to buffer particles during sorting for particles coming from tiles located on the same MPI rank, and another `nchunks_pad` chunks to buffer particles going to or coming from tiles located on different MPI ranks. The user can specify `nchunks_pad` directly, or can specify this as a fraction of the number of chunks per tile via `nchunks_pad_frac`. However, in practice, the default values should be appropriate for most applications. See further description of each below.

Related parameters are `ihole_size_frac` and `ihole_size`, which control the size of a buffer called "ihole" [1] used in the particle sorting algorithm. `ihole_size_frac` represents the number of particles which fit into ihole as a fraction of the total particle capacity of the chunk pool. For example, if `ihole_size_frac = 0.05`, then additional memory will be allocated to fit 5% of the particle capacity of the chunk pool; this memory is equally distributed among all tiles. The ihole itself does not use chunks, but because the required ihole size and `chunk_pool_size` both scale with the total number of particles, it is convenient to represent the size of the ihole as a fraction of the total chunk pool memory. The default value for `ihole_size_frac` of 0.05 is sufficient for most simulations. Alternatively, `ihole_size` can be specified. This parameter overrides `ihole_size_frac` to set the size of the ihole in absolute number of particles for each tile.

In the current implementation, if the ihole or chunks overflow during the particle sort, then the simulation crashes and prints an error informing the user to increase either `ihole_size_frac` or `nchunks_pad_frac` (or the absolute equivalents). In principle, overflows can be gracefully handled by using free chunks from the pool as additional memory within the particle sort. Normally, device memory cannot be reallocated with a CUDA kernel, but because the chunks are preallocated, they can be accessed and used arbitrarily. However, this form of dynamic error handling has not yet been implemented in OSIRIS.

In order to expedite communication for transfers between host and device, the code uses a "memory transfer manager" for batched data movements. This object buffers data and metadata for data movements in two fixed-size communication buffers allocated at initialization. One buffer resides on the host, the other, equal in size, resides on the device. The sending buffer is copied to the receiving buffer for movement between host and device whenever the sending buffer becomes full;  device-to-device copies between the device buffer and the original device memory addresses are performed via a kernel regardless of the movement direction. This construction greatly reduces the number of copies between host and device. The `com_buf_size` parameter controls the amount of memory allocated for the communication buffer in megabytes. The default value of 1GB should be sufficient for most applications.


## Using CUDA

To use CUDA, OSIRIS must be compiled CUDA support by setting `ENABLE_CUDA = 1` in the configuration script. See config/osiris_sys.perlmutter for an example of this and additional configuration options required.

For optimal performance:
* Thermal benchmarks indicate that at least $2^{18}$ total gridpoints should be used per GPU regardless of particles per cell.
* One should compile the code with the `nvcc` compiler flag `--maxrregcount`. On Perlmutter (NVIDIA A100 GPUs), best performance in 2D was achieved with `--maxrregcount 40`, and best performance in 3D was achieved with ``--maxrregcount 64``. This may vary by system/GPU.  See `config/osiris_sys.perlmutter.gpu` for an example.

### Input parameters

See `decks/cuda/` for example input decks.

* The cuda algorithm is currently implemented as a simulation mode which inherits from the `tiles` simulation mode. See `source/tiles/README.md` for a list of inherited input deck parameters. To activate the cuda algorithm, the following parameter must be added to the `simulation` section:
  * `algorithm=cuda,`

* Note on setting `node_number` in the `node_conf` section: The product of this number gives the total number of MPI ranks for the simulation. MPI ranks will be assigned to physical GPUs depending on how many are available. While it is possible to overload GPUs by setting the total number of MPI ranks to greater than the number of physical GPUs available, it is recommended that the user set the number of MPI ranks equal to the number of physical GPUs available (ie such that each MPI rank has its own GPU). This will give best performance. Also, the default particle buffer parameters above are set under this assumption. Running with more than one MPI rank per GPU can occasionally be useful for development. 

* Note on setting `tile_number` in the `node_conf` section: Maximum tile size is limited by the use of CUDA shared memory in the `dudt` kernel. Tiles must be small enough such that the E and B field arrays fit into shared memory. Therefore, maximum tile size is a function of simulation dimensionality, interpolation, precision (single or double), and the maximum shared memory per block on the GPU in question. Assuming tiles with sides of equal length (I.e., square-shaped in 2D, cube-shaped in 3D) the maximum tile size, in units of cells, per dimension is:

  ​        $n_x = (\frac{shmemsize}{2\*3\*precision})^{1/dims} - (2\*interp+1)$,

  where $shmemsize$ is the maximum shared memory size in bytes per block; $precision$ is 4 bytes for single precision, 8 bytes for double  precision; $dims$ is the simulation dimensionality, 1, 2 or 3; and $interp$ is 1 for linear, 2 for quadratic, 3 for cubic, and 4 for quartic interpolation. The factor of 2 in the denominator of the first term comes from the fact that there are 2 field arrays (E and B) that must fit into shared memory. The factor of 3 in the denominator of the first term comes from the fact that there are 3 field components. The second term is the number of guard cells.

  The maximum shared memory accessible per block ($shmemsize$) for CUDA GPUs of a given compute capability can be found in the [CUDA C programming guide](https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html?highlight=shared%20memory#compute-capability-8-x). The compute capability of different GPUs can be found at https://developer.nvidia.com/cuda-gpus. Alternatively, one can determine the max shmem per block through the terminal by compiling and running the `deviceQuery` program from the cuda samples https://github.com/NVIDIA/cuda-samples/tree/master/Samples/1_Utilities/deviceQuery.  

  For example, GPU nodes at Perlmutter at NERSC use NVIDIA A100 GPUs. NVIDIA A100 GPUs have compute capability 8.0. They have maximum shared memory per block of 163KB. 

* When using `algorithm=cuda` there is a different default value for `init_type` in the `species` namelist section: init_type, character(*), default = “transposed”. This is essentially the standard initialization, but tranposed such that particles which are close together in simulation space are scattered in memory, which gives better performance on the GPU by reducing memory collisions in the current deposit (`init_type="standard"` lays out particles in memory such that particles which are close in space are close in memory, which gives better performance on CPUs due to memory localization).

As described above, the default values of the rest of the parameters listed below should suffice, and are recommended, for most applications. 

* The following parameters have been added to the `node_conf` section:
  * `chunk_size`, integer, default = 512. Size of a chunk in units of number of particles. Must be greater than or equal to the block size, currently set to 512 (see `os-param-cuda.h`).
  * `chunk_pool_size`, integer, default = -1. Number of chunks allocated per MPI rank, in units of number of chunks. If not specified 75% of the available memory on the GPU will be allocated for chunks. This should be appropriate for most applications.
  * `com_buf_size`, integer, default = 1024. Size of buffer for CPU-GPU communication allocated for each MPI rank, in units of MB (so default is about 1GB). The default should be sufficient for most applications.
  * `ihole_size_frac`, real, default = 0.05. Amount of memory allocated for the ihole buffers for each tile specified as a fraction of the particle capacity of the chunk pool (`chunk_size`*`chunk_pool_size`). (Since this memory is distributed equally between all tiles, this parameter can also be thought of in terms of the approximate number of particles per tile.)
  * `ihole_size`, integer, default = -1. Optionally, one can specify the size of the ihole buffer for each tile in units of number of particles. If set, this overrides the value of `ihole_size_frac`.
  * `nchunks_pad_frac`, real, default = 0.05. Specifies, as a fraction of number of chunks per tile, (1) how many empty chunks each tile should have to buffer particles moving locally during the sort, and (2) how many chunks each tile should have to buffer particles moving externally (via MPI) during the sort.
  * `nchunks_pad`, integer, default = -1. Optionally, one can specify this in units of number of chunks. If set, this overrides the value of `nchunks_pad_frac`. 
  * `nblocks_requested`, integer, default = 10000. Number of GPU compute blocks requested. This default value should be appropriate for most applications.

## References

[1] Decyk, Viktor K., and Tajendra V. Singh. "Particle-in-cell algorithms for emerging computer architectures." *Computer Physics Communications* 185.3 (2014): 708-719./
