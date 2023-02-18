#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void expression_kernel(double *__restrict__ A, double const *__restrict__ w_0);
static void expression_kernel(double *__restrict__ A, double const *__restrict__ w_0)
{
  double const t0[7 * 3] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.0, 0.3333333333333334, 0.3333333333333333, 0.3333333333333333 };

  for (int32_t p0 = 0; p0 <= 6; ++p0)
    for (int32_t i = 0; i <= 2; ++i)
      A[p0] = A[p0] + t0[3 * p0 + i] * w_0[1 + 2 * i];
}

void wrap_expression_kernel(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ dat1, int32_t const *__restrict__ map0, int32_t const *__restrict__ map1)
{
  double t0[7];
  double t1[3 * 2];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    {
      int32_t const i9 = 0;

      for (int32_t i10 = 0; i10 <= 6; ++i10)
      {
        int32_t const i11 = 0;

        t0[i10] = 0.0;
      }
    }
    {
      int32_t const i12 = 0;

      for (int32_t i13 = 0; i13 <= 2; ++i13)
        for (int32_t i14 = 0; i14 <= 1; ++i14)
          t1[2 * i13 + i14] = dat1[2 * map1[3 * n + i13] + i14];
    }
    expression_kernel(&(t0[0]), &(t1[0]));
    {
      int32_t const i7 = 0;

      {
        int32_t const i8 = 0;

        for (int32_t i6 = 0; i6 <= 6; ++i6)
          dat0[map0[7 * n + i6]] = t0[i6];
      }
    }
  }
}