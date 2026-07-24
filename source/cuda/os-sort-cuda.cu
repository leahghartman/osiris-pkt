#include <stdio.h>
#include <inttypes.h>
#include <cuda_runtime.h>

#include "os-spec-push-cuda.h"
#include "os-param-cuda.h"


// --------------------------------------------------------------------------------------
// Perform exclusive parallel prefix sum on int array idata which is of size BLOCK_SIZE
// Note: Could implement avoiding bank conficts if this kernel is heavy
//       Tho it isn't. Still, we should run like it is godzilla. Ahhhhh
// --------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__device__ void prescan(int *const idata) {

  int tx = threadIdx.x;
  int offset = 1;

  // build sum in place up the tree
  for (int d = BLOCK_SIZE>>1; d > 0; d >>= 1) {
    __syncthreads();
    if (tx < d) {
      int ai = offset*(2*tx+1)-1;
      int bi = offset*(2*tx+2)-1;
      idata[bi] += idata[ai];
    }
    offset *= 2;
  }

  // clear the last element
  if (tx == 0) {
    idata[BLOCK_SIZE - 1] = 0;
  }

  // traverse down tree & build scan
  for (int d = 1; d < BLOCK_SIZE; d *= 2) {
    offset >>= 1;
    __syncthreads();
    if (tx < d) {
      int ai = offset*(2*tx+1)-1;
      int bi = offset*(2*tx+2)-1;
      int t = idata[ai];
      idata[ai] = idata[bi];
      idata[bi] += t;
    }
  }

}



// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__device__ int fill_holes( t_chunk **chunk, t_chunk **chunk_dest,
                           int num_par_proc_tot,
                           int part_offset, int n_holes_filled, t_f0<int> ihole,
                           int chunk_size, int nchunks_til, int *err ) {

  int tx = threadIdx.x;

  int n_holes_filled_temp;
  int chunk_idx;
  int part_idx_dest;

  // loop over chunks of particles
  int num_par_proc_chunk = min(num_par_proc_tot, chunk_size - part_offset);
  int part_idx = part_offset + tx;

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );
  t_f1<p_k_part> x_dest = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p_dest = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q_dest = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix_dest = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  while ( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    n_holes_filled_temp = n_holes_filled;

    // loop over particles in a chunk, BLOCK_SIZE at a time, move particles into holes
    while ( part_idx - part_offset < num_par_proc_chunk ) {

      // get index of hole
      part_idx_dest = ihole.fd( n_holes_filled_temp + tx );

      // get pointer to destination chunk buffers
      chunk_idx   = part_idx_dest / chunk_size;

      // crash if we try to access an invalid chunk
      if ( chunk_idx >= nchunks_til ) {
        // *err = p_err_nchunks_max;
        atomicExch( err, p_err_nchunks_max );
        return;
      }

      x_dest.g  = (p_k_part *) ( *(chunk_dest+chunk_idx) + x_offset );
      p_dest.g  = (p_k_part *) ( *(chunk_dest+chunk_idx) + p_offset );
      q_dest.g  = (p_k_part *) ( *(chunk_dest+chunk_idx) + q_offset );
      ix_dest.g = (int *)      ( *(chunk_dest+chunk_idx) + ix_offset);

      // write particle info to destination chunk

      for( int j=0; j<p_x_dim; j++ ) {
        x_dest.fd( part_idx_dest % chunk_size, j )  = x.fd( part_idx, j );
        ix_dest.fd( part_idx_dest % chunk_size, j ) = ix.fd( part_idx, j );
      }
      for( int j=0; j<p_p_dim; j++ ) {
        p_dest.fd( part_idx_dest % chunk_size, j ) = p.fd( part_idx, j );
      }
      q_dest.fd( part_idx_dest % chunk_size ) = q.fd( part_idx );

      // increment holes filled temp var
      n_holes_filled_temp += BLOCK_SIZE;

      part_idx += BLOCK_SIZE;

    } // end loop over chunk of particles

    // increment holes filled
    // this correctly handles case where we process fewer than BLOCK_SIZE number of threads
    // make it into the while loop on the last iteration for a given chunk
    // in that case n_holes_filled_temp overestimates how many holes we actually fileed
    n_holes_filled += num_par_proc_chunk;

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

    // now we will be starting at the beginning of a new chunk
    part_offset = 0;
    num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
    part_idx = tx;

  }

  return n_holes_filled;

}




// --------------------------------------------------------------------------------------
// Colonel Sortington: Sort particles by tile
// --------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void sort( t_chunk ***chunks, t_chunk ***chunks_ext,
                      int *nchunks_max, int *nchunks,
                      int *nchunks_max_ext, int *nchunks_ext,
                      int *my_nx_p_, int *block_start_tiles,
                      int *block_start_chunks, int *block_num_par, int chunk_size,
                      t_f0<int> *ihole_, int ihole_size, int *n_hole,
                      t_f1<int> *neighbor_til_id, t_f2<int> *shift, int *n_send_ext,
                      int *n_recv_loc, int *num_par, int *err ){

  // Block, thread index
  int bx = blockIdx.x;
  int tx = threadIdx.x;

  int my_til = block_start_tiles[bx];
  int chunk_idx = block_start_chunks[bx];
  int num_par_proc_tot = block_num_par[bx];

  t_chunk **chunk, **chunk_dest;

  chunk = chunks[my_til] + chunk_idx;
  int *my_nx_p = my_nx_p_ + (p_x_dim*my_til);
  t_f0<int> ihole = ihole_[my_til];

  int num_par_proc_chunk;
  int part_idx;
  int part_idx_dest;
  int chunk_idx_dest;

  // executable statements

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );
  t_f1<p_k_part> x_dest = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p_dest = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q_dest = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix_dest = t_f1<int>( chunk_size, p_x_dim );

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

      // deterimine direction particle is going
      int dir = 0;
      int fac = 1;
      for( int i=0; i<p_x_dim; i++ ) {
        dir += fac * ( ( ix.fd(part_idx,i) < 1 ) + 2 * ( ix.fd(part_idx,i) > my_nx_p[i] ) );
        fac *= 3;
      }

      // particle has left the bounds of this tile, move it
      if ( dir != 0 ) {

        // atomically update n_hole, write idx of moving particle to ihole
        int cur_hole = atomicAdd( n_hole + my_til, 1 );

        // check if we have enough space in ihole for this fresh new hole hot off the press
        if ( cur_hole > ihole_size ) {
          // *err = p_err_pierre_ihole;
          atomicExch( err, p_err_pierre_ihole );
          return;
        }

        // write particle index (0 indexed per tile) to the hole buffer
        ihole.fd( cur_hole ) = part_idx + chunk_idx*chunk_size;

        int dest_til = neighbor_til_id->fd( dir, my_til );

        // particle is going to tile on another node or crossing physical boundary
        if ( dest_til == -1 ) {

          // atomically update index for where we'll write the particle (0 indexed per tile)
          part_idx_dest = atomicAdd( n_send_ext + my_til, 1 );

          // get destination chunk
          chunk_idx_dest = part_idx_dest / chunk_size;

          // TODO do this
          // we've run (ran?) out of chunks and will need to get more from the pool...
          // for now just crash
          if ( chunk_idx_dest >= nchunks_ext[my_til] ) {
            // *err = p_err_nchunks_max_ext;
            atomicExch( err, p_err_nchunks_max_ext );
            return;

          // // we've run out of chunks and will need to grab more
          // } else if ( chunk_idx_dest >= nchunks_ext[my_til] ) {

          //   // for now just get as many chunks as we need (should just be 1, right?)
          //   nchunks_get = chunk_idx_dest + 1 - nchunks_ext[my_til];

          //   int ret = chunk_pool->get_chunks_device( nchunks_get, chunks_ext[my_til]
          //                                                        + nchunks_ext[my_til] );
          //   if ( ! ret ) {
          //     *err = p_err_nchunks_avail;
          //     return;
          //   }

          }

          chunk_dest = chunks_ext[my_til] + chunk_idx_dest;

        // particle going to another tile on this node
        } else {

          // atomically update index for where we'll write the particle
          part_idx_dest = atomicAdd( n_recv_loc + dest_til, 1 ) + num_par[dest_til];

          // get destination chunk
          chunk_idx_dest = part_idx_dest / chunk_size;

          // TODO do this, see above also
          // we've run out of chunks and will need to get more from the pool...using a lock...
          // for now just crash
          if ( chunk_idx_dest >= nchunks[dest_til] ) {
            // *err = p_err_nchunks_max;
            atomicExch( err, p_err_nchunks_max );
            return;
          }

          chunk_dest = chunks[dest_til] + chunk_idx_dest;

        }

        x_dest.g  = (p_k_part *) ( *chunk_dest + x_offset );
        p_dest.g  = (p_k_part *) ( *chunk_dest + p_offset );
        q_dest.g  = (p_k_part *) ( *chunk_dest + q_offset );
        ix_dest.g = (int *)      ( *chunk_dest + ix_offset);

        // write particle information to destination chunk
        for( int i=0; i<p_x_dim; i++ ) {
          x_dest.fd( part_idx_dest % chunk_size, i )  = x.fd( part_idx, i );
          ix_dest.fd( part_idx_dest % chunk_size, i ) = ix.fd( part_idx, i )
                                                        + shift->fd( i, dir, my_til );
        }
        for( int i=0; i<p_p_dim; i++ ) {
          p_dest.fd( part_idx_dest % chunk_size, i )  = p.fd( part_idx, i );
        }
        q_dest.fd( part_idx_dest % chunk_size )     = q.fd( part_idx );

        // tag holes in particle array for easy identification in compactification step
        ix.fd( part_idx, 0 ) = p_hole;

      }

      part_idx += BLOCK_SIZE;

    } // end loop over chunk of particles

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;
    chunk_idx += 1;

  } // end loop ``while ( num_par_proc_tot > 0 )``

}


// --------------------------------------------------------------------------------------
// Function stub to call kernel
// --------------------------------------------------------------------------------------
void sort_wrapper( t_chunk ***chunks, t_chunk ***chunks_ext,
                   int *nchunks_max, int *nchunks,
                   int *nchunks_max_ext, int *nchunks_ext,
                   int *my_nx_p, int *block_start_tiles, int *block_start_chunks,
                   int *block_num_par, int chunk_size, t_f0<int> *ihole, int ihole_size,
                   int *n_hole, t_f1<int> *neighbor_til_id, t_f2<int> *shift,
                   int *n_send_ext, int *n_recv_loc, int *num_par, int nblocks,
                   int *err ){

  sort<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__ >>>( chunks, chunks_ext, nchunks_max, nchunks,
                                        nchunks_max_ext, nchunks_ext,
                                        my_nx_p, block_start_tiles,
                                        block_start_chunks, block_num_par, chunk_size,
                                        ihole, ihole_size, n_hole, neighbor_til_id,
                                        shift, n_send_ext, n_recv_loc, num_par, err );

  chkLastErr("sort");
#ifdef __DEBUG__
  chkErr(cudaDeviceSynchronize());
#endif

}


// --------------------------------------------------------------------------------------
// Compactification station: Fill holes in particle array
// --------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void compactify( t_chunk ***chunks, t_chunk ***chunks_ext,
                            int *nchunks_max, int *nchunks,
                            int *nchunks_max_ext, int *nchunks_ext,
                            int chunk_size,
                            t_f0<int> *ihole_, int *n_hole_, int *num_par_,
                            int *n_recv_ext_, int *n_recv_loc_, int *err ){

  // Block, thread index
  int til = blockIdx.x;
  int tx = threadIdx.x;

  int n_hole = n_hole_[til];
  int n_recv_loc = n_recv_loc_[til];
  int num_par = num_par_[til];
  int n_recv_ext = n_recv_ext_[til];

  t_f0<int> ihole = ihole_[til];

  int num_par_new = num_par + n_recv_loc + n_recv_ext - n_hole;

  t_chunk **chunk, **chunk_dest;
  chunk_dest = chunks[til];

  int part_idx, part_idx_dest;
  int chunk_idx;
  int part_offset;
  int num_par_proc_tot, num_par_proc_chunk;
  int n_holes_filled;

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );
  t_f1<p_k_part> x_dest = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p_dest = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q_dest = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix_dest = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  // if we lost particles net, the sort is more complicated
  if ( num_par_new < num_par ) {
    
    __shared__ int mask[BLOCK_SIZE];
    int mask_last;

    // ------------
    // remove invalid holes from ihole ( holes with idx > num_par_new )
    // ------------
    int i = tx;
    int n_valid_hole = 0;
    int ihole_i;
    int n_hole_iters = (n_hole+BLOCK_SIZE-1)/BLOCK_SIZE;

    for ( int j=0; j<n_hole_iters; j++ ) {

      mask[tx] = 0;

      // store current hole to avoid race condition below
      if ( i < n_hole) {
        ihole_i = ihole.fd(i);

        // construct buffer indicating which threads encountered a valid hole
        if ( ihole_i < num_par_new ) {
          mask[tx] = 1;
        }
      }

      // save the last element of array since this is lost in exclusive prefix scan
      __syncthreads();
      mask_last = mask[BLOCK_SIZE-1];

      // perform prefix scan on mask
      // This gives us a mapping from threads with a valid hole to a new index in ihole
      // (ie stream compaction)
      prescan<BLOCK_SIZE>( mask );
      __syncthreads();

      if ( i < n_hole) {
        // write location of valid hole to ihole
        if ( ihole_i < num_par_new ) {
          ihole.fd( n_valid_hole + mask[tx] ) = ihole_i;
        }

        // increment valid hole counter
        n_valid_hole += mask[BLOCK_SIZE - 1] + mask_last;
      }

      i += BLOCK_SIZE;
      __syncthreads();

    }

    n_holes_filled = 0;

    // ------------
    // move particles with idx between num_par_new and num_par
    // (particles which remained on this til)
    // ------------
    chunk_idx = num_par_new / chunk_size;
    part_offset = num_par_new % chunk_size;
    num_par_proc_tot = num_par - num_par_new;

    chunk = chunks[til] + chunk_idx;

    // loop over chunks of particles
    num_par_proc_chunk = min(num_par_proc_tot, chunk_size - part_offset);
    part_idx = part_offset + tx;

    // loop over chunks of particles until all particles in tile have been processed
    while( num_par_proc_tot > 0 ) {

      x.g  = (p_k_part *) ( *chunk + x_offset );
      p.g  = (p_k_part *) ( *chunk + p_offset );
      q.g  = (p_k_part *) ( *chunk + q_offset );
      ix.g = (int *)      ( *chunk + ix_offset);

      int chunk_iters = (num_par_proc_chunk+BLOCK_SIZE-1)/BLOCK_SIZE;

      // loop over particles in a chunk, BLOCK_SIZE at a time
      for ( int j=0; j<chunk_iters; j++ ) {

        // every thread 0's out mask
        __syncthreads();
        mask[tx] = 0;

        if ( part_idx - part_offset < num_par_proc_chunk ) {

          // construct buffer indicating which threads encountered an actual particle (instead of a hole)
          if ( ix.fd( part_idx, 0 ) != p_hole ) {
            mask[tx] = 1;
          }

        }

        // save the last element of array since this is lost in exclusive prefix scan
        __syncthreads();
        mask_last = mask[BLOCK_SIZE-1];

        // perform prefix scan on mask
        // This gives us a mapping from threads with an actual particle to an index in ihole
        // (ie stream compaction)
        prescan<BLOCK_SIZE>( mask );
        __syncthreads();

        if ( part_idx - part_offset < num_par_proc_chunk ) {

          // move particles into holes
          if ( ix.fd( part_idx, 0 ) != p_hole ) {

            // get index of hole (0 indexed per tile)
            part_idx_dest = ihole.fd( n_holes_filled + mask[tx] );

            // get pointer to destination chunk buffers
            chunk_idx = part_idx_dest / chunk_size;

            // crash if we try to access an invalid chunk
            if ( chunk_idx >= nchunks[til] ) {
              // *err = p_err_nchunks_max;
              atomicExch( err, p_err_nchunks_max );
              return;
            }

            x_dest.g  = (p_k_part *) ( *(chunk_dest+chunk_idx) + x_offset );
            p_dest.g  = (p_k_part *) ( *(chunk_dest+chunk_idx) + p_offset );
            q_dest.g  = (p_k_part *) ( *(chunk_dest+chunk_idx) + q_offset );
            ix_dest.g = (int *)      ( *(chunk_dest+chunk_idx) + ix_offset);

            // write particle info to destination chunk
            for ( int j=0; j<p_x_dim; j++ ) {
              x_dest.fd( part_idx_dest % chunk_size, j ) = x.fd( part_idx, j );
              ix_dest.fd( part_idx_dest % chunk_size, j ) = ix.fd( part_idx, j );
            }

            for ( int j=0; j<p_p_dim; j++ ) {
              p_dest.fd( part_idx_dest % chunk_size, j ) = p.fd( part_idx, j );
            }

            q_dest.fd( part_idx_dest % chunk_size ) = q.fd( part_idx );

          }

        }

        // increment valid hole counter
        n_holes_filled += mask[BLOCK_SIZE - 1] + mask_last;

        part_idx += BLOCK_SIZE;

      } // end loop over chunk of particles

      num_par_proc_tot -= num_par_proc_chunk;
      chunk++;

      part_offset = 0;
      num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
      part_idx = tx;

    }

    // ------------
    // move particles with idx > num_par (particles we received locally)
    // ------------
    chunk_idx = num_par / chunk_size;
    part_offset = num_par % chunk_size;
    num_par_proc_tot = n_recv_loc;

    chunk = chunks[til] + chunk_idx;

    n_holes_filled = fill_holes<BLOCK_SIZE>( chunk, chunk_dest, num_par_proc_tot,
                                             part_offset, n_holes_filled, ihole,
                                             chunk_size, nchunks[til], err );
    if ( *err!=0 ) return;

    // ------------
    // move particles we received from MPI
    // ------------
    chunk_idx = 0;
    part_offset = 0;
    num_par_proc_tot = n_recv_ext;

    chunk = chunks_ext[til] + chunk_idx;

    n_holes_filled = fill_holes<BLOCK_SIZE>( chunk, chunk_dest, num_par_proc_tot,
                                             part_offset, n_holes_filled, ihole,
                                             chunk_size, nchunks[til], err );
    if ( *err!=0 ) return;

  // else if (num_par_new >= num_par)
  // If we gained or maintained number of particles net then
  // the sort is a bit more straight forward
  // Because: all holes are valid, and we don't have to faff about with sorting guys
  // which remained and are mixed in with holes
  } else {

    // Buffer is fully packed up to num_par_new with guys we received locally,
    // so we only need to fill holes.  And, these holes will be filled by local recv's
    // and external recv's, in an order that would surprise you
    if ( num_par + n_recv_loc >= num_par_new ) {

      n_holes_filled = 0;

      // ------------
      // move all particles we received locally with idx > num_par_new into holes
      // ------------
      chunk_idx = num_par_new / chunk_size;
      part_offset = num_par_new % chunk_size;
      num_par_proc_tot = num_par + n_recv_loc - num_par_new;

      chunk = chunks[til] + chunk_idx;

      n_holes_filled = fill_holes<BLOCK_SIZE>( chunk, chunk_dest, num_par_proc_tot,
                                               part_offset, n_holes_filled, ihole,
                                               chunk_size, nchunks[til], err );
      if ( *err!=0 ) return;

      // ------------
      // move particles we received from MPI
      // ------------

      chunk_idx = 0;
      part_offset = 0;
      num_par_proc_tot = n_recv_ext;

      chunk = chunks_ext[til] + chunk_idx;

      n_holes_filled = fill_holes<BLOCK_SIZE>( chunk, chunk_dest, num_par_proc_tot,
                                               part_offset, n_holes_filled, ihole,
                                               chunk_size, nchunks[til], err );
      if ( *err!=0 ) return;

    // Buffer is not fully packed up to num_par_new (ie we didn't recv that many
    // local guys, relative to how many guys we lost/how many we recvd externally).
    // In this case, we need to both fill holes, and spaces at the end of the buffer
    // only with the guys we recv'd from external comms
    } else {

      n_holes_filled = 0;

      // ------------
      // move particles we received from MPI into holes
      // ------------
      chunk_idx = 0;
      part_offset = 0;
      num_par_proc_tot = n_hole;

      chunk = chunks_ext[til] + chunk_idx;

      n_holes_filled = fill_holes<BLOCK_SIZE>( chunk, chunk_dest, num_par_proc_tot,
                                               part_offset, n_holes_filled, ihole,
                                               chunk_size, nchunks[til], err );
      if ( *err!=0 ) return;

      // ------------
      // move the rest of the particles we received from MPI into the end of the buffer
      // ------------
      int part_dest_offset = num_par + n_recv_loc;

      chunk_idx = n_hole / chunk_size;
      part_offset = n_hole % chunk_size;
      num_par_proc_tot = n_recv_ext - n_hole;

      chunk = chunks_ext[til] + chunk_idx;

      // loop over chunks of particles
      num_par_proc_chunk = min(num_par_proc_tot, chunk_size - part_offset);
      part_idx = part_offset + tx;

      while ( num_par_proc_tot > 0 ) {

        x.g  = (p_k_part *) ( *chunk + x_offset );
        p.g  = (p_k_part *) ( *chunk + p_offset );
        q.g  = (p_k_part *) ( *chunk + q_offset );
        ix.g = (int *)      ( *chunk + ix_offset);

        // get index of hole
        part_idx_dest = part_dest_offset + tx;

        // loop over particles in a chunk, BLOCK_SIZE at a time
        while ( part_idx - part_offset < num_par_proc_chunk ) {

          // move particles into holes

          // get pointer to destination chunk buffers
          chunk_idx = part_idx_dest / chunk_size;

          // crash if we try to access an invalid chunk
          if ( chunk_idx >= nchunks[til] ) {
            // *err = p_err_nchunks_max;
            atomicExch( err, p_err_nchunks_max );
            return;
          }
      
          x_dest.g  = (p_k_part *) ( *(chunk_dest+chunk_idx) + x_offset );
          p_dest.g  = (p_k_part *) ( *(chunk_dest+chunk_idx) + p_offset );
          q_dest.g  = (p_k_part *) ( *(chunk_dest+chunk_idx) + q_offset );
          ix_dest.g = (int *)      ( *(chunk_dest+chunk_idx) + ix_offset);

          // write particle info to destination chunk
          for( int j=0; j<p_x_dim; j++ ) {
            x_dest.fd( part_idx_dest % chunk_size, j ) = x.fd( part_idx, j );
            ix_dest.fd( part_idx_dest % chunk_size, j ) = ix.fd( part_idx, j );
          }
          for( int j=0; j<p_p_dim; j++ ) {
            p_dest.fd( part_idx_dest % chunk_size, j ) = p.fd( part_idx, j );
          }
          q_dest.fd( part_idx_dest % chunk_size ) = q.fd( part_idx );

          // increment destination particle index
          part_idx_dest += BLOCK_SIZE;

          part_idx += BLOCK_SIZE;

        } // end loop over chunk of particles

        part_dest_offset += num_par_proc_chunk;

        num_par_proc_tot -= num_par_proc_chunk;
        chunk++;

        part_offset = 0;
        num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
        part_idx = tx;

      }

    }

  }

  // update num particles on this tile
  if ( tx == 0 ) {
    num_par_[til] = num_par_new;
  }

}


// --------------------------------------------------------------------------------------
// Function stub to call kernel
// --------------------------------------------------------------------------------------
void compactify_wrapper( t_chunk ***chunks, t_chunk ***chunks_ext,
                         int *nchunks_max, int *nchunks,
                         int *nchunks_max_ext, int *nchunks_ext,
                         int chunk_size, t_f0<int> *ihole, int *n_hole, int *n_recv_ext,
                         int *num_par, int *n_recv_loc, int nblocks, int *err ){

  compactify<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__ >>>( chunks, chunks_ext,
                                              nchunks_max, nchunks,
                                              nchunks_max_ext, nchunks_ext,
                                              chunk_size, ihole,
                                              n_hole, num_par, n_recv_ext, n_recv_loc, err );

  chkLastErr("compactify");
#ifdef __DEBUG__
  chkErr(cudaDeviceSynchronize());
#endif

}
