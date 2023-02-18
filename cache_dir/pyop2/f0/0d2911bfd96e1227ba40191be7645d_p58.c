#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void form0_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2, double const *__restrict__ w_3, double const *__restrict__ w_4, double const *__restrict__ w_5, double const *__restrict__ w_6, double const *__restrict__ w_7, double const *__restrict__ w_8, double const *__restrict__ w_9);
static void form0_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2, double const *__restrict__ w_3, double const *__restrict__ w_4, double const *__restrict__ w_5, double const *__restrict__ w_6, double const *__restrict__ w_7, double const *__restrict__ w_8, double const *__restrict__ w_9)
{
  double t0;
  double t1;
  double t10;
  double t11;
  double t12;
  double const t13[7] = { 0.025, 0.025, 0.025, 0.06666666666666667, 0.06666666666666667, 0.06666666666666667, 0.225 };
  double const t14[7 * 7] = { -2.999999999999995, -1.0000000000000018, 0.0, 0.0, 5.24580540853904e-15, 3.9999999999999982, -5.181446385205573e-15, 1.0, 3.000000000000002, -1.1827932877306292e-15, 7.214931710316083e-15, 2.288682496747173e-15, -4.0, -1.0434330856368233e-14, 1.0000000000000018, -0.9999999999999946, -2.0016171622437674e-15, 4.000000000000007, -3.9999999999999907, -2.8690489112035944e-15, -1.8021493354097816e-14, 0.24999999999999853, 0.24999999999999986, -0.7499999999999998, 4.999999999999997, 1.0, 1.0, -6.749999999999997, -0.2500000000000002, -0.2500000000000012, 0.7500000000000006, -1.000000000000002, -5.000000000000002, -1.0, 6.750000000000004, -0.9999999999999988, 0.9999999999999978, 0.0, 0.0, -1.0770084712089784e-15, 0.0, 0.0, -0.33333333333333304, 0.3333333333333331, 0.0, 1.3333333333333326, -1.3333333333333324, 0.0, 0.0 };
  double const t15[7 * 7] = { -3.0000000000000013, 2.397800093331626e-15, -1.0000000000000018, -6.335379736374461e-15, 4.000000000000007, 7.800626541904116e-15, -8.881784197001252e-15, 0.9999999999999907, 4.745532949067312e-15, -1.0000000000000033, 4.0000000000000036, 2.2354338585325662e-14, -3.9999999999999987, -1.865174681370263e-14, 0.999999999999995, 6.123268552586174e-15, 3.0000000000000053, -7.545336504964141e-15, -3.9999999999999925, 0.0, -6.8832617168986165e-15, 0.2499999999999923, -0.7499999999999921, 0.2500000000000012, 4.999999999999992, 1.0000000000000164, 1.0, -6.750000000000011, -1.0000000000000053, 6.165839134505678e-15, 1.0000000000000024, -8.463143844880196e-15, 1.0604936603056867e-14, 3.990066140152413e-15, -8.8527487250726e-15, -0.25000000000000444, 0.7500000000000072, -0.2500000000000003, -1.0000000000000089, -0.9999999999999881, -5.0, 6.749999999999996, -0.33333333333333903, 6.698141318264943e-15, 0.3333333333333343, 1.3333333333333246, 1.1810485724339764e-14, -1.3333333333333321, -5.870663219447476e-15 };
  double const t16[7 * 7] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
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
  double t26;
  double t27;
  double t28;
  double t29;
  double t3;
  double t30;
  double t31;
  double t32;
  double t33;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;

  t0 = 1.0 / w_8[0];
  t1 = 1.0 / w_9[0];
  t2 = -0.9999999999999988 * coords[0] + 0.9999999999999988 * coords[2];
  t3 = -1.000000000000001 * coords[1] + coords[5];
  t4 = -1.000000000000001 * coords[0] + coords[4];
  t5 = -0.9999999999999988 * coords[1] + 0.9999999999999988 * coords[3];
  t6 = t2 * t3 + -1.0 * t4 * t5;
  t7 = 1.0 / t6;
  t8 = t2 * t7;
  t9 = t7 * -1.0 * t4;
  t10 = t7 * -1.0 * t5;
  t11 = t3 * t7;
  t12 = fabs(t6);
  for (int32_t ip = 0; ip <= 6; ++ip)
  {
    t26 = 0.0;
    t25 = 0.0;
    t24 = 0.0;
    t23 = 0.0;
    t22 = 0.0;
    t21 = 0.0;
    t20 = 0.0;
    t19 = 0.0;
    t18 = 0.0;
    t17 = 0.0;
    for (int32_t i = 0; i <= 6; ++i)
    {
      t17 = t17 + t16[7 * ip + i] * w_3[1 + 2 * i];
      t18 = t18 + t16[7 * ip + i] * w_3[2 * i];
    }
    for (int32_t i_0 = 0; i_0 <= 6; ++i_0)
    {
      t19 = t19 + t16[7 * ip + i_0] * w_4[i_0];
      t20 = t20 + t16[7 * ip + i_0] * w_2[i_0];
      t21 = t21 + t16[7 * ip + i_0] * w_6[i_0];
      t22 = t22 + t16[7 * ip + i_0] * w_5[i_0];
      t23 = t23 + t16[7 * ip + i_0] * w_1[i_0];
      t24 = t24 + t15[7 * ip + i_0] * w_2[i_0];
      t25 = t25 + t14[7 * ip + i_0] * w_2[i_0];
      t26 = t26 + t16[7 * ip + i_0] * w_0[i_0];
    }
    t27 = t13[ip] * t12;
    t28 = t26 * t26;
    t29 = t11 * t25 + t10 * t24;
    t30 = t9 * t25 + t8 * t24;
    t31 = t27 * (-1.0 * t28 * (t29 * t10 + t30 * t8) + -1.0 * (t18 * t10 + t17 * t8));
    t32 = t27 * (-1.0 * t28 * (t29 * t11 + t30 * t9) + -1.0 * (t18 * t11 + t17 * t9));
    t33 = t27 * (-1.0 * (-1.0 * 2.0 * t20 + t19) * t0 + -1.0 * (t22 + t21) * t1 * -1.0 * t19 + t28 * t23 * w_7[0] + -1.0 * t20 * t22 * t21);
    for (int32_t j = 0; j <= 6; ++j)
      A[j] = A[j] + t15[7 * ip + j] * t31 + t16[7 * ip + j] * t33 + t14[7 * ip + j] * t32;
  }
}

void wrap_form0_cell_integral_otherwise(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ dat1, double const *__restrict__ dat2, double const *__restrict__ dat3, double const *__restrict__ dat4, double const *__restrict__ dat5, double const *__restrict__ dat6, double const *__restrict__ dat7, double const *__restrict__ dat8, double const *__restrict__ glob0, double const *__restrict__ glob1, double const *__restrict__ glob2, int32_t const *__restrict__ map0, int32_t const *__restrict__ map1)
{
  double t0[7];
  double t1[3 * 2];
  double t2[7];
  double t3[7];
  double t4[7];
  double t5[7 * 2];
  double t6[7];
  double t7[7];
  double t8[7];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    for (int32_t i40 = 0; i40 <= 6; ++i40)
    {
      {
        int32_t const i52 = 0;

        {
          int32_t const i53 = 0;

          t8[i40] = dat8[map0[7 * n + i40]];
        }
      }
      {
        int32_t const i50 = 0;

        {
          int32_t const i51 = 0;

          t7[i40] = dat7[map0[7 * n + i40]];
        }
      }
      {
        int32_t const i48 = 0;

        {
          int32_t const i49 = 0;

          t6[i40] = dat6[map0[7 * n + i40]];
        }
      }
      {
        int32_t const i46 = 0;

        for (int32_t i47 = 0; i47 <= 1; ++i47)
          t5[2 * i40 + i47] = dat5[2 * map0[7 * n + i40] + i47];
      }
      {
        int32_t const i44 = 0;

        {
          int32_t const i45 = 0;

          t4[i40] = dat4[map0[7 * n + i40]];
        }
      }
      {
        int32_t const i42 = 0;

        {
          int32_t const i43 = 0;

          t3[i40] = dat3[map0[7 * n + i40]];
        }
      }
      {
        int32_t const i39 = 0;

        {
          int32_t const i41 = 0;

          t2[i40] = dat2[map0[7 * n + i40]];
        }
      }
    }
    {
      int32_t const i36 = 0;

      for (int32_t i37 = 0; i37 <= 2; ++i37)
        for (int32_t i38 = 0; i38 <= 1; ++i38)
          t1[2 * i37 + i38] = dat1[2 * map1[3 * n + i37] + i38];
    }
    {
      int32_t const i33 = 0;

      for (int32_t i34 = 0; i34 <= 6; ++i34)
      {
        int32_t const i35 = 0;

        t0[i34] = 0.0;
      }
    }
    form0_cell_integral_otherwise(&(t0[0]), &(t1[0]), &(t2[0]), &(t3[0]), &(t4[0]), &(t5[0]), &(t6[0]), &(t7[0]), &(t8[0]), &(glob0[0]), &(glob1[0]), &(glob2[0]));
    {
      int32_t const i31 = 0;

      {
        int32_t const i32 = 0;

        for (int32_t i30 = 0; i30 <= 6; ++i30)
          dat0[map0[7 * n + i30]] = dat0[map0[7 * n + i30]] + t0[i30];
      }
    }
  }
}