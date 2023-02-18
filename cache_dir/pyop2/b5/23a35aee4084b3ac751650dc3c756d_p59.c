#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void expression_kernel(double *__restrict__ A, double const *__restrict__ w_0);
static void expression_kernel(double *__restrict__ A, double const *__restrict__ w_0)
{
  for (int32_t p0 = 0; p0 <= 6; ++p0)
    for (int32_t p1 = 0; p1 <= 1; ++p1)
      A[2 * p0 + p1] = A[2 * p0 + p1] + w_0[p1];
}

void wrap_expression_kernel(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ glob0, int32_t const *__restrict__ map0)
{
  double t0[7 * 2];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    {
      int32_t const i7 = 0;

      for (int32_t i8 = 0; i8 <= 6; ++i8)
        for (int32_t i9 = 0; i9 <= 1; ++i9)
          t0[2 * i8 + i9] = 0.0;
    }
    expression_kernel(&(t0[0]), &(glob0[0]));
    for (int32_t i5 = 0; i5 <= 1; ++i5)
    {
      int32_t const i6 = 0;

      for (int32_t i4 = 0; i4 <= 6; ++i4)
        dat0[2 * map0[7 * n + i4] + i5] = t0[2 * i4 + i5];
    }
  }
}