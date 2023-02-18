#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void expression_kernel(double *__restrict__ A, double const *__restrict__ w_0);
static void expression_kernel(double *__restrict__ A, double const *__restrict__ w_0)
{
  double const t0[7 * 7] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };

  for (int32_t p0 = 0; p0 <= 6; ++p0)
    for (int32_t i = 0; i <= 6; ++i)
      A[p0] = A[p0] + t0[7 * p0 + i] * w_0[i];
}

void wrap_expression_kernel(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ dat1, int32_t const *__restrict__ map0)
{
  double t0[7];
  double t1[7];

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

      for (int32_t i13 = 0; i13 <= 6; ++i13)
      {
        int32_t const i14 = 0;

        t1[i13] = dat1[map0[7 * n + i13]];
      }
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