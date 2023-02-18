#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void expression_kernel(double *__restrict__ A, double const *__restrict__ w_0);
static void expression_kernel(double *__restrict__ A, double const *__restrict__ w_0)
{
  double const t0[6 * 7] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 };

  for (int32_t p0 = 0; p0 <= 5; ++p0)
    for (int32_t p1 = 0; p1 <= 1; ++p1)
      for (int32_t i = 0; i <= 6; ++i)
        A[2 * p0 + p1] = A[2 * p0 + p1] + t0[7 * p0 + i] * w_0[p1 + 2 * i];
}

void wrap_expression_kernel(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ dat1, int32_t const *__restrict__ map0, int32_t const *__restrict__ map1)
{
  double t0[6 * 2];
  double t1[7 * 2];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    {
      int32_t const i9 = 0;

      for (int32_t i10 = 0; i10 <= 5; ++i10)
        for (int32_t i11 = 0; i11 <= 1; ++i11)
          t0[2 * i10 + i11] = 0.0;
    }
    {
      int32_t const i12 = 0;

      for (int32_t i13 = 0; i13 <= 6; ++i13)
        for (int32_t i14 = 0; i14 <= 1; ++i14)
          t1[2 * i13 + i14] = dat1[2 * map1[7 * n + i13] + i14];
    }
    expression_kernel(&(t0[0]), &(t1[0]));
    for (int32_t i7 = 0; i7 <= 1; ++i7)
    {
      int32_t const i8 = 0;

      for (int32_t i6 = 0; i6 <= 5; ++i6)
        dat0[2 * map0[6 * n + i6] + i7] = t0[2 * i6 + i7];
    }
  }
}