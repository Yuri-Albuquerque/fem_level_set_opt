#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void form0_interior_facet_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2, uint32_t const *__restrict__ facet);
static void form0_interior_facet_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2, uint32_t const *__restrict__ facet)
{
  double const t0[3 * 2] = { -1.0, 1.0, 0.0, 1.0, 1.0, 0.0 };
  double t1;
  double const t10[3 * 2] = { 1.0, 1.0, -1.0, -0.0, -0.0, -1.0 };
  double t11;
  double t12;
  double t13;
  double t14;
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t2;
  double t20;
  double t21;
  double t22;
  double t23;
  double t24;
  double t25;
  double const t26[3 * 2 * 7] = { 0.0, 0.4553418012614795, -0.12200846792814622, 0.6666666666666667, 0.0, 0.0, 0.0, 0.0, -0.12200846792814618, 0.4553418012614795, 0.6666666666666667, 0.0, 0.0, 0.0, 0.4553418012614795, 0.0, -0.12200846792814622, 0.0, 0.6666666666666667, 0.0, 0.0, -0.12200846792814621, 0.0, 0.4553418012614795, 0.0, 0.6666666666666667, 0.0, 0.0, 0.4553418012614795, -0.12200846792814618, 0.0, 0.0, 0.0, 0.6666666666666667, 0.0, -0.12200846792814622, 0.4553418012614795, 0.0, 0.0, 0.0, 0.6666666666666667, 0.0 };
  double t27;
  double t28;
  double t29;
  double t3;
  double t30;
  double t31;
  double t32;
  double t33;
  double t34;
  double t35;
  double const t36[1] = { 1.0 };
  double t4;
  double t5;
  double t6;
  double t7;
  double const t8[2] = { 0.5, 0.5 };
  double t9;

  t33 = 0.0;
  t1 = -0.9999999999999988 * coords[0] + 0.9999999999999988 * coords[2];
  t2 = -1.0 * coords[0] + coords[4];
  t3 = t0[2 * facet[0]] * t1 + t0[2 * facet[0] + 1] * t2;
  t4 = -0.9999999999999988 * coords[1] + 0.9999999999999988 * coords[3];
  t5 = -1.0 * coords[1] + coords[5];
  t6 = t0[2 * facet[0]] * t4 + t0[2 * facet[0] + 1] * t5;
  t7 = sqrt(t3 * t3 + t6 * t6);
  t9 = 1.0 / (t1 * t5 + -1.0 * t2 * t4);
  t11 = t10[2 * facet[0]] * t5 * t9 + t10[2 * facet[0] + 1] * t9 * -1.0 * t4;
  t12 = t10[2 * facet[0]] * t9 * -1.0 * t2 + t10[2 * facet[0] + 1] * t1 * t9;
  t13 = 1.0 / sqrt(t11 * t11 + t12 * t12);
  t14 = t12 * t13;
  t15 = t11 * t13;
  t16 = -0.9999999999999988 * coords[6] + 0.9999999999999988 * coords[8];
  t17 = -1.0 * coords[7] + coords[11];
  t18 = -1.0 * coords[6] + coords[10];
  t19 = -0.9999999999999988 * coords[7] + 0.9999999999999988 * coords[9];
  t20 = 1.0 / (t16 * t17 + -1.0 * t18 * t19);
  t21 = t10[2 * facet[1]] * t17 * t20 + t10[2 * facet[1] + 1] * t20 * -1.0 * t19;
  t22 = t10[2 * facet[1]] * t20 * -1.0 * t18 + t10[2 * facet[1] + 1] * t16 * t20;
  t23 = 1.0 / sqrt(t21 * t21 + t22 * t22);
  t24 = t22 * t23;
  t25 = t21 * t23;
  for (int32_t ip = 0; ip <= 1; ++ip)
  {
    t30 = 0.0;
    t29 = 0.0;
    t28 = 0.0;
    t27 = 0.0;
    for (int32_t i = 0; i <= 6; ++i)
    {
      t27 = t27 + t26[14 * facet[0] + 7 * ip + i] * w_0[1 + 2 * i];
      t28 = t28 + t26[14 * facet[0] + 7 * ip + i] * w_0[2 * i];
      t29 = t29 + t26[14 * facet[1] + 7 * ip + i] * w_0[15 + 2 * i];
      t30 = t30 + t26[14 * facet[1] + 7 * ip + i] * w_0[14 + 2 * i];
    }
    t31 = t25 * 1.2 * t30 + t24 * 1.2 * t29;
    t32 = t15 * 1.2 * t28 + t14 * 1.2 * t27;
    t33 = t33 + (-1.0 * w_2[1] * 0.5 * (fabs(t31) + t31) + w_2[0] * 0.5 * (fabs(t32) + t32)) * t8[ip] * t7;
  }
  t34 = t33 * -1.0 * w_1[0];
  t35 = w_1[0] * t33;
  {
    int32_t const j = 0;

    A[0] = A[0] + t36[0] * t35;
    A[1] = A[1] + t36[0] * t34;
  }
}

void wrap_form0_interior_facet_integral_otherwise(int32_t const start, int32_t const end, double *__restrict__ dat1, double const *__restrict__ dat2, double const *__restrict__ dat3, double const *__restrict__ glob0, double const *__restrict__ dat4, uint32_t const *__restrict__ dat0, int32_t const *__restrict__ map0, int32_t const *__restrict__ map1, int32_t const *__restrict__ map2)
{
  double t0[2];
  double t1[6 * 2];
  double t2[14 * 2];
  double t3[2];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    {
      int32_t const i26 = 0;

      for (int32_t i27 = 0; i27 <= 1; ++i27)
      {
        int32_t const i28 = 0;

        t3[i27] = dat4[map0[2 * n + i27]];
      }
    }
    {
      int32_t const i23 = 0;

      for (int32_t i24 = 0; i24 <= 13; ++i24)
        for (int32_t i25 = 0; i25 <= 1; ++i25)
          t2[2 * i24 + i25] = dat3[2 * map2[14 * n + i24] + i25];
    }
    {
      int32_t const i20 = 0;

      for (int32_t i21 = 0; i21 <= 5; ++i21)
        for (int32_t i22 = 0; i22 <= 1; ++i22)
          t1[2 * i21 + i22] = dat2[2 * map1[6 * n + i21] + i22];
    }
    {
      int32_t const i17 = 0;

      for (int32_t i18 = 0; i18 <= 1; ++i18)
      {
        int32_t const i19 = 0;

        t0[i18] = 0.0;
      }
    }
    form0_interior_facet_integral_otherwise(&(t0[0]), &(t1[0]), &(t2[0]), &(glob0[0]), &(t3[0]), &(dat0[2 * n]));
    {
      int32_t const i15 = 0;

      {
        int32_t const i16 = 0;

        for (int32_t i14 = 0; i14 <= 1; ++i14)
          dat1[map0[2 * n + i14]] = dat1[map0[2 * n + i14]] + t0[i14];
      }
    }
  }
}