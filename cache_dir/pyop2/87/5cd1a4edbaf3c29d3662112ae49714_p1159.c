#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void form0_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2, double const *__restrict__ w_3);
static void form0_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2, double const *__restrict__ w_3)
{
  double t0;
  double t1;
  double t10;
  double t11;
  double const t12[3] = { 0.16666666666666666, 0.16666666666666666, 0.16666666666666666 };
  double const t13[3 * 7] = { -1.4166666666666663, -0.08333333333333492, 0.2500000000000001, -0.33333333333333415, -1.6666666666666665, 1.0000000000000013, 2.250000000000001, 0.3333333333333329, -0.3333333333333326, 0.0, 2.6666666666666647, -2.6666666666666647, 0.0, 0.0, 0.08333333333333481, 1.4166666666666667, -0.25000000000000017, 1.666666666666667, 0.33333333333333415, -1.0000000000000016, -2.2500000000000013 };
  double const t14[3 * 7] = { -1.416666666666672, 0.2500000000000055, -0.083333333333333, -0.3333333333333406, 1.0000000000000109, -1.6666666666666645, 2.249999999999994, 0.08333333333332746, -0.2499999999999929, 1.4166666666666685, 1.6666666666666565, -0.9999999999999883, 0.33333333333333637, -2.2500000000000075, 0.33333333333332726, 7.175496405181828e-15, -0.33333333333333304, 2.666666666666659, 1.3986955577605752e-14, -2.6666666666666665, -8.15791481085019e-15 };
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t2;
  double const t20[1] = { 1.0 };
  double t3;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;

  t19 = 0.0;
  t0 = w_2[0] * -1.0 * w_1[0];
  t1 = -0.9999999999999988 * coords[0] + 0.9999999999999988 * coords[2];
  t2 = -1.0 * coords[1] + coords[5];
  t3 = -1.0 * coords[0] + coords[4];
  t4 = -0.9999999999999988 * coords[1] + 0.9999999999999988 * coords[3];
  t5 = t1 * t2 + -1.0 * t3 * t4;
  t6 = 1.0 / t5;
  t7 = t1 * t6;
  t8 = t6 * -1.0 * t3;
  t9 = t6 * -1.0 * t4;
  t10 = t2 * t6;
  t11 = fabs(t5);
  for (int32_t ip = 0; ip <= 2; ++ip)
  {
    t18 = 0.0;
    t17 = 0.0;
    t16 = 0.0;
    t15 = 0.0;
    for (int32_t i = 0; i <= 6; ++i)
    {
      t15 = t15 + t14[7 * ip + i] * w_0[1 + 2 * i];
      t16 = t16 + t13[7 * ip + i] * w_0[1 + 2 * i];
      t17 = t17 + t14[7 * ip + i] * w_0[2 * i];
      t18 = t18 + t13[7 * ip + i] * w_0[2 * i];
    }
    t19 = t19 + (t0 * (0.5120000000000001 * (t10 * t18 + t9 * t17) + 0.5120000000000001 * (t8 * t16 + t7 * t15)) + w_3[0]) * t12[ip] * t11;
  }
  {
    int32_t const j = 0;

    A[0] = A[0] + t20[0] * t19;
  }
}

void wrap_form0_cell_integral_otherwise(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ dat1, double const *__restrict__ dat2, double const *__restrict__ dat3, double const *__restrict__ glob0, double const *__restrict__ dat4, int32_t const *__restrict__ map0, int32_t const *__restrict__ map1, int32_t const *__restrict__ map2)
{
  double t0[1];
  double t1[3 * 2];
  double t2[7 * 2];
  double t3[1];
  double t4[1];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    {
      int32_t const i31 = 0;

      {
        int32_t const i32 = 0;

        {
          int32_t const i33 = 0;

          t4[0] = dat4[map0[n]];
        }
      }
    }
    {
      int32_t const i28 = 0;

      {
        int32_t const i29 = 0;

        {
          int32_t const i30 = 0;

          t3[0] = dat3[map0[n]];
        }
      }
    }
    {
      int32_t const i25 = 0;

      for (int32_t i26 = 0; i26 <= 6; ++i26)
        for (int32_t i27 = 0; i27 <= 1; ++i27)
          t2[2 * i26 + i27] = dat2[2 * map2[7 * n + i26] + i27];
    }
    {
      int32_t const i22 = 0;

      for (int32_t i23 = 0; i23 <= 2; ++i23)
        for (int32_t i24 = 0; i24 <= 1; ++i24)
          t1[2 * i23 + i24] = dat1[2 * map1[3 * n + i23] + i24];
    }
    {
      int32_t const i19 = 0;

      {
        int32_t const i20 = 0;

        {
          int32_t const i21 = 0;

          t0[0] = 0.0;
        }
      }
    }
    form0_cell_integral_otherwise(&(t0[0]), &(t1[0]), &(t2[0]), &(t3[0]), &(glob0[0]), &(t4[0]));
    {
      int32_t const i17 = 0;

      {
        int32_t const i18 = 0;

        {
          int32_t const i16 = 0;

          dat0[map0[n]] = dat0[map0[n]] + t0[0];
        }
      }
    }
  }
}