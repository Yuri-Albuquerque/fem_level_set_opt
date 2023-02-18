#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void form0_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2, double const *__restrict__ w_3);
static void form0_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2, double const *__restrict__ w_3)
{
  double t0;
  double t1;
  double t2;
  double const t3[7] = { 0.025, 0.025, 0.025, 0.06666666666666667, 0.06666666666666667, 0.06666666666666667, 0.225 };
  double const t4[7 * 7] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
  double t5;
  double t6;
  double t7;

  t0 = 1.0 / w_2[0];
  t1 = 1.0 / w_3[0];
  t2 = fabs((-0.9999999999999988 * coords[0] + 0.9999999999999988 * coords[2]) * (-1.000000000000001 * coords[1] + coords[5]) + -1.0 * (-1.000000000000001 * coords[0] + coords[4]) * (-0.9999999999999988 * coords[1] + 0.9999999999999988 * coords[3]));
  for (int32_t ip = 0; ip <= 6; ++ip)
  {
    t6 = 0.0;
    t5 = 0.0;
    for (int32_t i = 0; i <= 6; ++i)
    {
      t5 = t5 + t4[7 * ip + i] * w_1[i];
      t6 = t6 + t4[7 * ip + i] * w_0[i];
    }
    t7 = ((t6 + t5) * t1 + t0) * t3[ip] * t2;
    for (int32_t k = 0; k <= 6; ++k)
      A[k] = A[k] + t7 * t4[7 * ip + k] * t4[7 * ip + k];
  }
}

void wrap_form0_cell_integral_otherwise(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ dat1, double const *__restrict__ dat2, double const *__restrict__ dat3, double const *__restrict__ glob0, double const *__restrict__ glob1, int32_t const *__restrict__ map0, int32_t const *__restrict__ map1)
{
  double t0[7];
  double t1[3 * 2];
  double t2[7];
  double t3[7];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    for (int32_t i24 = 0; i24 <= 6; ++i24)
    {
      {
        int32_t const i26 = 0;

        {
          int32_t const i27 = 0;

          t3[i24] = dat3[map0[7 * n + i24]];
        }
      }
      {
        int32_t const i23 = 0;

        {
          int32_t const i25 = 0;

          t2[i24] = dat2[map0[7 * n + i24]];
        }
      }
    }
    {
      int32_t const i20 = 0;

      for (int32_t i21 = 0; i21 <= 2; ++i21)
        for (int32_t i22 = 0; i22 <= 1; ++i22)
          t1[2 * i21 + i22] = dat1[2 * map1[3 * n + i21] + i22];
    }
    {
      int32_t const i17 = 0;

      for (int32_t i18 = 0; i18 <= 6; ++i18)
      {
        int32_t const i19 = 0;

        t0[i18] = 0.0;
      }
    }
    form0_cell_integral_otherwise(&(t0[0]), &(t1[0]), &(t2[0]), &(t3[0]), &(glob0[0]), &(glob1[0]));
    {
      int32_t const i15 = 0;

      {
        int32_t const i16 = 0;

        for (int32_t i14 = 0; i14 <= 6; ++i14)
          dat0[map0[7 * n + i14]] = dat0[map0[7 * n + i14]] + t0[i14];
      }
    }
  }
}