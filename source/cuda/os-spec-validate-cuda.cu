#include <stdio.h>
#include <inttypes.h>
#include <cuda_runtime.h>

#include "os-spec-push-cuda.h"


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void validate( t_chunk ***chunks, int *nchunks_, int *my_nx_p_, int *block_start_tiles,
                          int *block_start_chunks, int *block_num_par, int chunk_size,
                          bool over_, int move_num0, int move_num1, int move_num2, int *err,
                          int my_aid ){

  // Block, thread index
  int bx = blockIdx.x;
  int tx = threadIdx.x;

  int my_til = block_start_tiles[bx];
  int chunk_idx = block_start_chunks[bx];
  int num_par_proc_tot = block_num_par[bx];

  t_chunk **chunk;

  chunk = chunks[my_til] + chunk_idx;
  int nchunks = nchunks_[my_til];
  int *my_nx_p = my_nx_p_ + (p_x_dim*my_til);

  int num_par_proc_chunk;
  int part_idx;

  int move_num[3] = { move_num0, move_num1, move_num2 };

  // sometimes we allow for 1 cell overflow (e.g. before update boundary)
  int over = 0;
  if (over_==true) over = 1;

  // get parameters to handle chunks
  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  // loop over chunks of particles until all particles in tile have been processed
  while( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
    part_idx = tx;

    // loop over particles in a chunk, BLOCK_SIZE at a time
    while ( part_idx < num_par_proc_chunk ) {

      // validate positions
      for( int i=0; i<p_x_dim; i++ ) {

        int ilb = 1 - move_num[i] - over;
        int iub = my_nx_p[i] + over;

        if ( x.fd(part_idx,i) < p_k_part(-0.5) or x.fd(part_idx,i) >= p_k_part(0.5) or
             ix.fd(part_idx,i) < ilb or ix.fd(part_idx,i) > iub ) {
          printf("[Node %d], tile %d: Bad particle position! At index %d of %d, on chunk %d of %d. "
                 "x[%d]=%f, ix[%d]=%d\n",
                 my_aid, my_til, part_idx, chunk_size, chunk_idx, nchunks,
                 i, x.fd(part_idx,i), i, ix.fd(part_idx,i));
          *err = 1;
        }

      }

      // validate momenta
      for( int i=0; i<p_p_dim; i++ ) {
        if ( !isfinite(p.fd(part_idx,i)) ) {
          printf("[Node %d], tile %d: Bad particle momentum! At index %d of %d, on chunk %d of %d. "
                 "p[%d]=%f \n",
                 my_aid, my_til, part_idx, chunk_size, chunk_idx, nchunks,
                 i, p.fd(part_idx,i) );
          *err = 1;
        }
      }

      // validate charge
      if ( !isfinite(q.fd(part_idx)) ) {
        printf("[Node %d], tile %d: Bad particle charge! At index %d of %d, on chunk %d of %d. "
               "q=%f\n",
               my_aid, my_til, part_idx, chunk_size, chunk_idx, nchunks,
               q.fd(part_idx) );
        *err = 1;
      }

      part_idx += BLOCK_SIZE;

    } // end loop over chunk of particles ``while ( part_idx < num_par_proc_chunk )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;
    chunk_idx += 1;

  } // end loop ``while ( num_par_proc_tot > 0 )``

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void validate_wrapper( t_chunk ***chunks, int *nchunks, int *my_nx_p_, int *block_start_tiles,
                       int *block_start_chunks, int *block_num_par, int chunk_size,
                       bool over_, int move_num0, int move_num1, int move_num2,
                       int *err, int my_aid, int nblocks){

  validate<__BLOCK_SIZE__><<< nblocks, __BLOCK_SIZE__ >>>( chunks, nchunks, my_nx_p_, block_start_tiles,
                                           block_start_chunks, block_num_par, chunk_size,
                                           over_, move_num0, move_num1, move_num2, err,
                                           my_aid);

  chkLastErr("validate_wrapper");
#ifdef __DEBUG__
  chkErr(cudaDeviceSynchronize());
#endif

}
