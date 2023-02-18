#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void form0_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2);
static void form0_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2)
{
  double t0;
  double t1;
  double const t2[7] = { 0.025, 0.025, 0.025, 0.06666666666666667, 0.06666666666666667, 0.06666666666666667, 0.225 };
  double const t3[7 * 7] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
  double t4;
  double t5;
  double t6;

  t0 = 1.0 / w_2[0];
  t1 = fabs((-0.9999999999999988 * coords[0] + 0.9999999999999988 * coords[2]) * (-1.000000000000001 * coords[1] + coords[5]) + -1.0 * (-1.000000000000001 * coords[0] + coords[4]) * (-0.9999999999999988 * coords[1] + 0.9999999999999988 * coords[3]));
  for (int32_t ip = 0; ip <= 6; ++ip)
  {
    t5 = 0.0;
    t4 = 0.0;
    for (int32_t i = 0; i <= 6; ++i)
    {
      t4 = t4 + t3[7 * ip + i] * w_1[i];
      t5 = t5 + t3[7 * ip + i] * w_0[i];
    }
    t6 = ((t5 + t4) * 10000.0 + t0) * t2[ip] * t1;
    for (int32_t k = 0; k <= 6; ++k)
      A[k] = A[k] + t6 * t3[7 * ip + k] * t3[7 * ip + k];
  }
}

void wrap_form0_cell_integral_otherwise(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ dat1, double const *__restrict__ dat2, double const *__restrict__ dat3, double const *__restrict__ glob0, int32_t const *__restrict__ map0, int32_t const *__restrict__ map1)
{
  double t0[7];
  double t1[3 * 2];
  double t2[7];
  double t3[7];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    for (int32_t i23 = 0; i23 <= 6; ++i23)
    {
      {
        int32_t const i25 = 0;

        {
          int32_t const i26 = 0;

          t3[i23] = dat3[map0[7 * n + i23]];
        }
      }
      {
        int32_t const i22 = 0;

        {
          int32_t const i24 = 0;

          t2[i23] = dat2[map0[7 * n + i23]];
        }
      }
    }
    {
      int32_t const i19 = 0;

      for (int32_t i20 = 0; i20 <= 2; ++i20)
        for (int32_t i21 = 0; i21 <= 1; ++i21)
          t1[2 * i20 + i21] = dat1[2 * map1[3 * n + i20] + i21];
    }
    {
      int32_t const i16 = 0;

      for (int32_t i17 = 0; i17 <= 6; ++i17)
      {
        int32_t const i18 = 0;

        t0[i17] = 0.0;
      }
    }
    form0_cell_integral_otherwise(&(t0[0]), &(t1[0]), &(t2[0]), &(t3[0]), &(glob0[0]));
    {
      int32_t const i14 = 0;

      {
        int32_t const i15 = 0;

        for (int32_t i13 = 0; i13 <= 6; ++i13)
          dat0[map0[7 * n + i13]] = dat0[map0[7 * n + i13]] + t0[i13];
      }
    }
  }
}