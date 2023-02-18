#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void expression_kernel(double *__restrict__ A, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2);
static void expression_kernel(double *__restrict__ A, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2)
{
  double t0;
  double t1;
  double t2;
  double const t3[7 * 3] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.0, 0.3333333333333334, 0.3333333333333333, 0.3333333333333333 };
  double t4;
  double t5;

  t0 = -1.0 * w_2[0];
  t1 = -1.0 * w_1[1];
  t2 = -1.0 * w_1[0];
  for (int32_t p0 = 0; p0 <= 6; ++p0)
  {
    t5 = 0.0;
    t4 = 0.0;
    for (int32_t i = 0; i <= 2; ++i)
    {
      t4 = t4 + t3[3 * p0 + i] * w_0[1 + 2 * i];
      t5 = t5 + t3[3 * p0 + i] * w_0[2 * i];
    }
    A[p0] = A[p0] + exp((pow(t5 + t2, 2.0) + pow(t4 + t1, 2.0)) * t0);
  }
}

void wrap_expression_kernel(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ dat1, double const *__restrict__ glob0, double const *__restrict__ glob1, int32_t const *__restrict__ map0, int32_t const *__restrict__ map1)
{
  double t0[7];
  double t1[3 * 2];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    {
      int32_t const i14 = 0;

      for (int32_t i15 = 0; i15 <= 2; ++i15)
        for (int32_t i16 = 0; i16 <= 1; ++i16)
          t1[2 * i15 + i16] = dat1[2 * map1[3 * n + i15] + i16];
    }
    {
      int32_t const i11 = 0;

      for (int32_t i12 = 0; i12 <= 6; ++i12)
      {
        int32_t const i13 = 0;

        t0[i12] = 0.0;
      }
    }
    expression_kernel(&(t0[0]), &(t1[0]), &(glob0[0]), &(glob1[0]));
    {
      int32_t const i9 = 0;

      for (int32_t i8 = 0; i8 <= 6; ++i8)
      {
        int32_t const i10 = 0;

        dat0[map0[7 * n + i8]] = t0[i8];
      }
    }
  }
}