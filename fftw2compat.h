#ifndef FFTW2COMPAT_H
#define FFTW2COMPAT_H

#include <stdlib.h>
#include <fftw3-mpi.h>
#include <mpi.h>

#ifndef FFTW_REAL_TO_COMPLEX
#define FFTW_REAL_TO_COMPLEX 0
#endif
#ifndef FFTW_COMPLEX_TO_REAL
#define FFTW_COMPLEX_TO_REAL 1
#endif
#ifndef FFTW_NORMAL_ORDER
#define FFTW_NORMAL_ORDER 0
#endif

typedef double fftw_real;

typedef struct {
  fftw_plan plan;
  int direction;
} rfftwnd_mpi_plan_data;

typedef rfftwnd_mpi_plan_data *rfftwnd_mpi_plan;

static inline rfftwnd_mpi_plan rfftw3d_mpi_create_plan(MPI_Comm comm,
                                                       int nx,
                                                       int ny,
                                                       int nz,
                                                       int direction,
                                                       unsigned flags)
{
  ptrdiff_t local_n0, local_0_start;
  ptrdiff_t alloc_local;
  fftw_real *buffer;
  fftw_complex *cbuffer;
  rfftwnd_mpi_plan wrapper;

  alloc_local = fftw_mpi_local_size_3d(nx, ny, nz / 2 + 1, comm,
                                       &local_n0, &local_0_start);
  buffer = fftw_alloc_real(2 * alloc_local);
  if(buffer == NULL)
    return NULL;
  cbuffer = (fftw_complex *) buffer;

  wrapper = malloc(sizeof(*wrapper));
  if(wrapper == NULL)
    {
      fftw_free(buffer);
      return NULL;
    }

  wrapper->direction = direction;
  if(direction == FFTW_REAL_TO_COMPLEX)
    wrapper->plan = fftw_mpi_plan_dft_r2c_3d(nx, ny, nz, buffer, cbuffer, comm, flags);
  else
    wrapper->plan = fftw_mpi_plan_dft_c2r_3d(nx, ny, nz, cbuffer, buffer, comm, flags);

  fftw_free(buffer);

  if(wrapper->plan == NULL)
    {
      free(wrapper);
      return NULL;
    }

  return wrapper;
}

static inline void rfftwnd_mpi_local_sizes(rfftwnd_mpi_plan plan,
                                           int *local_nx,
                                           int *local_x_start,
                                           int *local_ny_after_transpose,
                                           int *local_y_start_after_transpose,
                                           int *total_size,
                                           int nx,
                                           int ny,
                                           int nz,
                                           MPI_Comm comm)
{
  ptrdiff_t local_n0, local_0_start;
  ptrdiff_t alloc_local;
  (void)plan;
  alloc_local = fftw_mpi_local_size_3d(nx, ny, nz / 2 + 1, comm,
                                       &local_n0, &local_0_start);
  *local_nx = (int)local_n0;
  *local_x_start = (int)local_0_start;
  if(local_ny_after_transpose)
    *local_ny_after_transpose = 0;
  if(local_y_start_after_transpose)
    *local_y_start_after_transpose = 0;
  *total_size = (int)(2 * alloc_local);
}

static inline void rfftwnd_mpi(rfftwnd_mpi_plan plan,
                               int howmany,
                               fftw_real *data,
                               fftw_real *workspace,
                               int order)
{
  (void)howmany;
  (void)workspace;
  (void)order;
  if(plan->direction == FFTW_REAL_TO_COMPLEX)
    fftw_mpi_execute_dft_r2c(plan->plan, data, (fftw_complex *) data);
  else
    fftw_mpi_execute_dft_c2r(plan->plan, (fftw_complex *) data, data);
}

static inline void rfftwnd_mpi_destroy_plan(rfftwnd_mpi_plan plan)
{
  if(plan != NULL)
    {
      if(plan->plan != NULL)
        fftw_destroy_plan(plan->plan);
      free(plan);
    }
}

#endif
