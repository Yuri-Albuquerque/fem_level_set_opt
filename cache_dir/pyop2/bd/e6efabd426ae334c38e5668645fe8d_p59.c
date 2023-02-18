#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void expression_kernel(double *__restrict__ A, double const *__restrict__ w_0, double const *__restrict__ w_1);
static void expression_kernel(double *__restrict__ A, double const *__restrict__ w_0, double const *__restrict__ w_1)
{
  double const t0[7 * 7] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
  double t1;
  double t2;

  for (int32_t p0 = 0; p0 <= 6; ++p0)
  {
    t2 = 0.0;
    t1 = 0.0;
    for (int32_t i = 0; i <= 6; ++i)
    {
      t1 = t1 + t0[7 * p0 + i] * w_1[i];
      t2 = t2 + t0[7 * p0 + i] * w_0[i];
    }
    A[p0] = A[p0] + t2 + t1;
  }
}

void wrap_expression_kernel(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ dat1, double const *__restrict__ dat2, int32_t const *__restrict__ map0)
{
  double t0[7];
  double t1[7];
  double t2[7];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    for (int32_t i16 = 0; i16 <= 6; ++i16)
    {
      {
        int32_t const i18 = 0;

        {
          int32_t const i19 = 0;

          t2[i16] = dat2[map0[7 * n + i16]];
        }
      }
      {
        int32_t const i15 = 0;

        {
          int32_t const i17 = 0;

          t1[i16] = dat1[map0[7 * n + i16]];
        }
      }
    }
    {
      int32_t const i12 = 0;

      for (int32_t i13 = 0; i13 <= 6; ++i13)
      {
        int32_t const i14 = 0;

        t0[i13] = 0.0;
      }
    }
    expression_kernel(&(t0[0]), &(t1[0]), &(t2[0]));
    {
      int32_t const i10 = 0;

      for (int32_t i9 = 0; i9 <= 6; ++i9)
      {
        int32_t const i11 = 0;

        dat0[map0[7 * n + i9]] = t0[i9];
      }
    }
  }
}