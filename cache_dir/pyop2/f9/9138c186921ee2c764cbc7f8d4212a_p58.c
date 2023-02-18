#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void par_loop_kernel(double const *__restrict__ f, double *__restrict__ f_max, double *__restrict__ f_min);
static void par_loop_kernel(double const *__restrict__ f, double *__restrict__ f_max, double *__restrict__ f_min)
{
  for (int32_t i = 0; i <= 2; ++i)
    for (int32_t d = 0; d <= 1; ++d)
    {
      f_min[d] = fmin(f_min[d], f[d + 2 * i]);
      f_max[d] = fmax(f_max[d], f[d + 2 * i]);
    }
}

void wrap_par_loop_kernel(int32_t const start, int32_t const end, double const *__restrict__ dat2, double *__restrict__ dat0, double *__restrict__ dat1, int32_t const *__restrict__ map1, int32_t const *__restrict__ map0)
{
  double t0[3 * 2];
  double t1[2];
  double t2[2];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    for (int32_t i19 = 0; i19 <= 1; ++i19)
    {
      {
        int32_t const i20 = 0;

        {
          int32_t const i21 = 0;

          t2[i19] = dat1[2 * map0[n] + i19];
        }
      }
      {
        int32_t const i17 = 0;

        {
          int32_t const i18 = 0;

          t1[i19] = dat0[2 * map0[n] + i19];
        }
      }
    }
    {
      int32_t const i14 = 0;

      for (int32_t i15 = 0; i15 <= 2; ++i15)
        for (int32_t i16 = 0; i16 <= 1; ++i16)
          t0[2 * i15 + i16] = dat2[2 * map1[3 * n + i15] + i16];
    }
    par_loop_kernel(&(t0[0]), &(t1[0]), &(t2[0]));
    for (int32_t i10 = 0; i10 <= 1; ++i10)
    {
      {
        int32_t const i9 = 0;

        {
          int32_t const i11 = 0;

          dat0[2 * map0[n] + i10] = fmax(dat0[2 * map0[n] + i10], t1[i10]);
        }
      }
      {
        int32_t const i13 = 0;

        {
          int32_t const i12 = 0;

          dat1[2 * map0[n] + i10] = fmin(dat1[2 * map0[n] + i10], t2[i10]);
        }
      }
    }
  }
}