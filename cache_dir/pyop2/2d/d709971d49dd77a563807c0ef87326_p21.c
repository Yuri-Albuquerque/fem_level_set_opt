#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void expression_kernel(double *__restrict__ A, double const *__restrict__ w_0);
static void expression_kernel(double *__restrict__ A, double const *__restrict__ w_0)
{
  for (int32_t p0 = 0; p0 <= 6; ++p0)
    A[p0] = A[p0] + w_0[0];
}

void wrap_expression_kernel(int32_t const start, int32_t const end, int32_t const *__restrict__ subset_indices, double *__restrict__ dat0, double const *__restrict__ glob0, int32_t const *__restrict__ map0)
{
  int32_t t0;
  double t1[7];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    t0 = subset_indices[n];
    {
      int32_t const i7 = 0;

      for (int32_t i8 = 0; i8 <= 6; ++i8)
      {
        int32_t const i9 = 0;

        t1[i8] = 0.0;
      }
    }
    expression_kernel(&(t1[0]), &(glob0[0]));
    {
      int32_t const i5 = 0;

      {
        int32_t const i6 = 0;

        for (int32_t i4 = 0; i4 <= 6; ++i4)
          dat0[map0[7 * t0 + i4]] = t1[i4];
      }
    }
  }
}