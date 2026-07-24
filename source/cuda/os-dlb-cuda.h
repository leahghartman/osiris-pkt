#ifndef OS_DLB_CUDA_H
#define OS_DLB_CUDA_H

#include "os-current-cuda.h"
#include "os-emf-cuda.h"

void copy_tils_stnry_wrapper( int *tils_stnry, int *tils_stnry_old, int n_tils_stnry,
                              t_current_c *jay, t_current_c *jay_temp,
                              t_emf_c *emf, t_emf_c *emf_temp, int nblocks );

#endif
