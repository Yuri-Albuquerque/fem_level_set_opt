#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void form0_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2, double const *__restrict__ w_3, double const *__restrict__ w_4);
static void form0_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2, double const *__restrict__ w_3, double const *__restrict__ w_4)
{
  double const t0[7 * 7] = { -3.0000000000000013, 2.397800093331626e-15, -1.0000000000000018, -6.335379736374461e-15, 4.000000000000007, 7.800626541904116e-15, -8.881784197001252e-15, 0.9999999999999907, 4.745532949067312e-15, -1.0000000000000033, 4.0000000000000036, 2.2354338585325662e-14, -3.9999999999999987, -1.865174681370263e-14, 0.999999999999995, 6.123268552586174e-15, 3.0000000000000053, -7.545336504964141e-15, -3.9999999999999925, 0.0, -6.8832617168986165e-15, 0.2499999999999923, -0.7499999999999921, 0.2500000000000012, 4.999999999999992, 1.0000000000000164, 1.0, -6.750000000000011, -1.0000000000000053, 6.165839134505678e-15, 1.0000000000000024, -8.463143844880196e-15, 1.0604936603056867e-14, 3.990066140152413e-15, -8.8527487250726e-15, -0.25000000000000444, 0.7500000000000072, -0.2500000000000003, -1.0000000000000089, -0.9999999999999881, -5.0, 6.749999999999996, -0.33333333333333903, 6.698141318264943e-15, 0.3333333333333343, 1.3333333333333246, 1.1810485724339764e-14, -1.3333333333333321, -5.870663219447476e-15 };
  double const t1[7 * 7] = { -2.999999999999995, -1.0000000000000018, 0.0, 0.0, 5.24580540853904e-15, 3.9999999999999982, -5.181446385205573e-15, 1.0, 3.000000000000002, -1.1827932877306292e-15, 7.214931710316083e-15, 2.288682496747173e-15, -4.0, -1.0434330856368233e-14, 1.0000000000000018, -0.9999999999999946, -2.0016171622437674e-15, 4.000000000000007, -3.9999999999999907, -2.8690489112035944e-15, -1.8021493354097816e-14, 0.24999999999999853, 0.24999999999999986, -0.7499999999999998, 4.999999999999997, 1.0, 1.0, -6.749999999999997, -0.2500000000000002, -0.2500000000000012, 0.7500000000000006, -1.000000000000002, -5.000000000000002, -1.0, 6.750000000000004, -0.9999999999999988, 0.9999999999999978, 0.0, 0.0, -1.0770084712089784e-15, 0.0, 0.0, -0.33333333333333304, 0.3333333333333331, 0.0, 1.3333333333333326, -1.3333333333333324, 0.0, 0.0 };
  double t10;
  double t11;
  double t12;
  double const t13[7] = { 0.025, 0.025, 0.025, 0.06666666666666667, 0.06666666666666667, 0.06666666666666667, 0.225 };
  double const t14[7 * 7] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
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
  double t3;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;

  t2 = -1.000000000000001 * coords[0] + coords[4];
  t3 = -0.9999999999999988 * coords[0] + 0.9999999999999988 * coords[2];
  t4 = -1.000000000000001 * coords[1] + coords[5];
  t5 = -0.9999999999999988 * coords[1] + 0.9999999999999988 * coords[3];
  t6 = t3 * t4 + -1.0 * t2 * t5;
  t7 = 1.0 / t6;
  t8 = t7 * -1.0 * t2;
  t9 = t4 * t7;
  t10 = t3 * t7;
  t11 = t7 * -1.0 * t5;
  t12 = fabs(t6);
  for (int32_t ip = 0; ip <= 6; ++ip)
  {
    t19 = 0.0;
    t18 = 0.0;
    t17 = 0.0;
    t16 = 0.0;
    t15 = 0.0;
    for (int32_t i = 0; i <= 6; ++i)
    {
      t15 = t15 + t14[7 * ip + i] * w_1[i];
      t16 = t16 + t14[7 * ip + i] * w_4[i];
      t17 = t17 + t14[7 * ip + i] * w_3[i];
      t18 = t18 + t14[7 * ip + i] * w_2[i];
      t19 = t19 + t14[7 * ip + i] * w_0[i];
    }
    t20 = t13[ip] * t12;
    t21 = -1.0 * pow(t19, 2.0);
    t22 = t20 * (t21 * (t18 * t11 + t17 * t10) + t10 * t16);
    t23 = t20 * (t21 * (t18 * t9 + t17 * t8) + t8 * t16);
    t24 = t20 * (t21 * (t18 * t10 + t15 * t11) + t11 * t16);
    t25 = t20 * (t21 * (t18 * t8 + t15 * t9) + t9 * t16);
    for (int32_t j0 = 0; j0 <= 6; ++j0)
    {
      A[2 * j0] = A[2 * j0] + t1[7 * ip + j0] * t25 + t0[7 * ip + j0] * t24;
      A[1 + 2 * j0] = A[1 + 2 * j0] + t1[7 * ip + j0] * t23 + t0[7 * ip + j0] * t22;
    }
  }
}

void wrap_form0_cell_integral_otherwise(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ dat1, double const *__restrict__ dat2, double const *__restrict__ dat3, double const *__restrict__ dat4, double const *__restrict__ dat5, double const *__restrict__ dat6, int32_t const *__restrict__ map0, int32_t const *__restrict__ map1)
{
  double t0[7 * 2];
  double t1[3 * 2];
  double t2[7];
  double t3[7];
  double t4[7];
  double t5[7];
  double t6[7];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    for (int32_t i31 = 0; i31 <= 6; ++i31)
    {
      {
        int32_t const i39 = 0;

        {
          int32_t const i40 = 0;

          t6[i31] = dat6[map0[7 * n + i31]];
        }
      }
      {
        int32_t const i37 = 0;

        {
          int32_t const i38 = 0;

          t5[i31] = dat5[map0[7 * n + i31]];
        }
      }
      {
        int32_t const i35 = 0;

        {
          int32_t const i36 = 0;

          t4[i31] = dat4[map0[7 * n + i31]];
        }
      }
      {
        int32_t const i33 = 0;

        {
          int32_t const i34 = 0;

          t3[i31] = dat3[map0[7 * n + i31]];
        }
      }
      {
        int32_t const i30 = 0;

        {
          int32_t const i32 = 0;

          t2[i31] = dat2[map0[7 * n + i31]];
        }
      }
    }
    {
      int32_t const i27 = 0;

      for (int32_t i28 = 0; i28 <= 2; ++i28)
        for (int32_t i29 = 0; i29 <= 1; ++i29)
          t1[2 * i28 + i29] = dat1[2 * map1[3 * n + i28] + i29];
    }
    {
      int32_t const i24 = 0;

      for (int32_t i25 = 0; i25 <= 6; ++i25)
        for (int32_t i26 = 0; i26 <= 1; ++i26)
          t0[2 * i25 + i26] = 0.0;
    }
    form0_cell_integral_otherwise(&(t0[0]), &(t1[0]), &(t2[0]), &(t3[0]), &(t4[0]), &(t5[0]), &(t6[0]));
    for (int32_t i22 = 0; i22 <= 1; ++i22)
    {
      int32_t const i23 = 0;

      for (int32_t i21 = 0; i21 <= 6; ++i21)
        dat0[2 * map0[7 * n + i21] + i22] = dat0[2 * map0[7 * n + i21] + i22] + t0[2 * i21 + i22];
    }
  }
}