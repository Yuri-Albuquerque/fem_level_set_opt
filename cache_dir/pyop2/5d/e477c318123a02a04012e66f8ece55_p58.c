#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void form0_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1);
static void form0_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1)
{
  double t0;
  double const t1[7] = { 0.025, 0.025, 0.025, 0.06666666666666667, 0.06666666666666667, 0.06666666666666667, 0.225 };
  double const t2[7 * 7] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
  double t3;
  double t4;
  double t5;

  t0 = fabs((-0.9999999999999988 * coords[0] + 0.9999999999999988 * coords[2]) * (-1.000000000000001 * coords[1] + coords[5]) + -1.0 * (-1.000000000000001 * coords[0] + coords[4]) * (-0.9999999999999988 * coords[1] + 0.9999999999999988 * coords[3]));
  for (int32_t ip = 0; ip <= 6; ++ip)
  {
    t4 = 0.0;
    t3 = 0.0;
    for (int32_t i = 0; i <= 6; ++i)
    {
      t3 = t3 + t2[7 * ip + i] * w_1[i];
      t4 = t4 + t2[7 * ip + i] * w_0[i];
    }
    t5 = t4 * t3 * t1[ip] * t0;
    for (int32_t j = 0; j <= 6; ++j)
      A[j] = A[j] + t2[7 * ip + j] * t5;
  }
}

void wrap_form0_cell_integral_otherwise(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ dat1, double const *__restrict__ dat2, double const *__restrict__ dat3, int32_t const *__restrict__ map0, int32_t const *__restrict__ map1)
{
  double t0[7];
  double t1[3 * 2];
  double t2[7];
  double t3[7];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    for (int32_t i22 = 0; i22 <= 6; ++i22)
    {
      {
        int32_t const i24 = 0;

        {
          int32_t const i25 = 0;

          t3[i22] = dat3[map0[7 * n + i22]];
        }
      }
      {
        int32_t const i21 = 0;

        {
          int32_t const i23 = 0;

          t2[i22] = dat2[map0[7 * n + i22]];
        }
      }
    }
    {
      int32_t const i18 = 0;

      for (int32_t i19 = 0; i19 <= 2; ++i19)
        for (int32_t i20 = 0; i20 <= 1; ++i20)
          t1[2 * i19 + i20] = dat1[2 * map1[3 * n + i19] + i20];
    }
    {
      int32_t const i15 = 0;

      for (int32_t i16 = 0; i16 <= 6; ++i16)
      {
        int32_t const i17 = 0;

        t0[i16] = 0.0;
      }
    }
    form0_cell_integral_otherwise(&(t0[0]), &(t1[0]), &(t2[0]), &(t3[0]));
    {
      int32_t const i13 = 0;

      {
        int32_t const i14 = 0;

        for (int32_t i12 = 0; i12 <= 6; ++i12)
          dat0[map0[7 * n + i12]] = dat0[map0[7 * n + i12]] + t0[i12];
      }
    }
  }
}