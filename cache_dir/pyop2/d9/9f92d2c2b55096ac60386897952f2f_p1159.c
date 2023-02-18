#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void form0_exterior_facet_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2, double const *__restrict__ w_3, uint32_t const *__restrict__ facet);
static void form0_exterior_facet_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2, double const *__restrict__ w_3, uint32_t const *__restrict__ facet)
{
  double const t0[3 * 2] = { 1.0, 1.0, -1.0, -0.0, -0.0, -1.0 };
  double t1;
  double t10;
  double const t11[3 * 2] = { -1.0, 1.0, 0.0, 1.0, 1.0, 0.0 };
  double t12;
  double t13;
  double t14;
  double const t15[2] = { 0.5, 0.5 };
  double const t16[3 * 2 * 7] = { 0.0, 0.4553418012614795, -0.12200846792814622, 0.6666666666666667, 0.0, 0.0, 0.0, 0.0, -0.12200846792814618, 0.4553418012614795, 0.6666666666666667, 0.0, 0.0, 0.0, 0.4553418012614795, 0.0, -0.12200846792814622, 0.0, 0.6666666666666667, 0.0, 0.0, -0.12200846792814621, 0.0, 0.4553418012614795, 0.0, 0.6666666666666667, 0.0, 0.0, 0.4553418012614795, -0.12200846792814618, 0.0, 0.0, 0.0, 0.6666666666666667, 0.0, -0.12200846792814622, 0.4553418012614795, 0.0, 0.0, 0.0, 0.6666666666666667, 0.0 };
  double t17;
  double t18;
  double t19;
  double t2;
  double t20;
  double const t21[1] = { 1.0 };
  double t3;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;

  t20 = 0.0;
  t1 = -0.9999999999999988 * coords[0] + 0.9999999999999988 * coords[2];
  t2 = -1.0 * coords[1] + coords[5];
  t3 = -1.0 * coords[0] + coords[4];
  t4 = -0.9999999999999988 * coords[1] + 0.9999999999999988 * coords[3];
  t5 = 1.0 / (t1 * t2 + -1.0 * t3 * t4);
  t6 = t0[2 * facet[0]] * t2 * t5 + t0[2 * facet[0] + 1] * t5 * -1.0 * t4;
  t7 = t0[2 * facet[0]] * t5 * -1.0 * t3 + t0[2 * facet[0] + 1] * t1 * t5;
  t8 = 1.0 / sqrt(t6 * t6 + t7 * t7);
  t9 = t7 * t8;
  t10 = t6 * t8;
  t12 = t11[2 * facet[0]] * t1 + t11[2 * facet[0] + 1] * t3;
  t13 = t11[2 * facet[0]] * t4 + t11[2 * facet[0] + 1] * t2;
  t14 = sqrt(t12 * t12 + t13 * t13);
  for (int32_t ip = 0; ip <= 1; ++ip)
  {
    t18 = 0.0;
    t17 = 0.0;
    for (int32_t i = 0; i <= 6; ++i)
    {
      t17 = t17 + t16[14 * facet[0] + 7 * ip + i] * w_0[1 + 2 * i];
      t18 = t18 + t16[14 * facet[0] + 7 * ip + i] * w_0[2 * i];
    }
    t19 = t10 * 0.6400000000000001 * t18 + t9 * 0.6400000000000001 * t17;
    t20 = t20 + (w_1[0] * (t19 < 0.0 ? t19 * w_2[0] : 0.0) + w_1[0] * (t19 > 0.0 ? t19 * w_3[0] : 0.0)) * t15[ip] * t14;
  }
  {
    int32_t const j = 0;

    A[0] = A[0] + t21[0] * t20;
  }
}

void wrap_form0_exterior_facet_integral_otherwise(int32_t const start, int32_t const end, double *__restrict__ dat1, double const *__restrict__ dat2, double const *__restrict__ dat3, double const *__restrict__ glob0, double const *__restrict__ glob1, double const *__restrict__ dat4, uint32_t const *__restrict__ dat0, int32_t const *__restrict__ map0, int32_t const *__restrict__ map1, int32_t const *__restrict__ map2)
{
  double t0[1];
  double t1[3 * 2];
  double t2[7 * 2];
  double t3[1];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    {
      int32_t const i27 = 0;

      {
        int32_t const i28 = 0;

        {
          int32_t const i29 = 0;

          t3[0] = dat4[map0[n]];
        }
      }
    }
    {
      int32_t const i24 = 0;

      for (int32_t i25 = 0; i25 <= 6; ++i25)
        for (int32_t i26 = 0; i26 <= 1; ++i26)
          t2[2 * i25 + i26] = dat3[2 * map2[7 * n + i25] + i26];
    }
    {
      int32_t const i21 = 0;

      for (int32_t i22 = 0; i22 <= 2; ++i22)
        for (int32_t i23 = 0; i23 <= 1; ++i23)
          t1[2 * i22 + i23] = dat2[2 * map1[3 * n + i22] + i23];
    }
    {
      int32_t const i18 = 0;

      {
        int32_t const i19 = 0;

        {
          int32_t const i20 = 0;

          t0[0] = 0.0;
        }
      }
    }
    form0_exterior_facet_integral_otherwise(&(t0[0]), &(t1[0]), &(t2[0]), &(glob0[0]), &(glob1[0]), &(t3[0]), &(dat0[n]));
    {
      int32_t const i16 = 0;

      {
        int32_t const i17 = 0;

        {
          int32_t const i15 = 0;

          dat1[map0[n]] = dat1[map0[n]] + t0[0];
        }
      }
    }
  }
}