#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void form0_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1);
static void form0_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1)
{
  double const t0[7 * 7] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
  double t1;
  double t10;
  double t11;
  double t12;
  double const t13[7 * 7] = { -2.999999999999995, -1.0000000000000018, 0.0, 0.0, 5.24580540853904e-15, 3.9999999999999982, -5.181446385205573e-15, 1.0, 3.000000000000002, -1.1827932877306292e-15, 7.214931710316083e-15, 2.288682496747173e-15, -4.0, -1.0434330856368233e-14, 1.0000000000000018, -0.9999999999999946, -2.0016171622437674e-15, 4.000000000000007, -3.9999999999999907, -2.8690489112035944e-15, -1.8021493354097816e-14, 0.24999999999999853, 0.24999999999999986, -0.7499999999999998, 4.999999999999997, 1.0, 1.0, -6.749999999999997, -0.2500000000000002, -0.2500000000000012, 0.7500000000000006, -1.000000000000002, -5.000000000000002, -1.0, 6.750000000000004, -0.9999999999999988, 0.9999999999999978, 0.0, 0.0, -1.0770084712089784e-15, 0.0, 0.0, -0.33333333333333304, 0.3333333333333331, 0.0, 1.3333333333333326, -1.3333333333333324, 0.0, 0.0 };
  double const t14[7 * 7] = { -3.0000000000000013, 2.397800093331626e-15, -1.0000000000000018, -6.335379736374461e-15, 4.000000000000007, 7.800626541904116e-15, -8.881784197001252e-15, 0.9999999999999907, 4.745532949067312e-15, -1.0000000000000033, 4.0000000000000036, 2.2354338585325662e-14, -3.9999999999999987, -1.865174681370263e-14, 0.999999999999995, 6.123268552586174e-15, 3.0000000000000053, -7.545336504964141e-15, -3.9999999999999925, 0.0, -6.8832617168986165e-15, 0.2499999999999923, -0.7499999999999921, 0.2500000000000012, 4.999999999999992, 1.0000000000000164, 1.0, -6.750000000000011, -1.0000000000000053, 6.165839134505678e-15, 1.0000000000000024, -8.463143844880196e-15, 1.0604936603056867e-14, 3.990066140152413e-15, -8.8527487250726e-15, -0.25000000000000444, 0.7500000000000072, -0.2500000000000003, -1.0000000000000089, -0.9999999999999881, -5.0, 6.749999999999996, -0.33333333333333903, 6.698141318264943e-15, 0.3333333333333343, 1.3333333333333246, 1.1810485724339764e-14, -1.3333333333333321, -5.870663219447476e-15 };
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t2;
  double t20;
  double t21;
  double t3;
  double t4;
  double t5;
  double t6;
  double const t7[7] = { 0.025, 0.025, 0.025, 0.06666666666666667, 0.06666666666666667, 0.06666666666666667, 0.225 };
  double t8;
  double t9;

  t1 = -0.9999999999999988 * coords[0] + 0.9999999999999988 * coords[2];
  t2 = -1.000000000000001 * coords[1] + coords[5];
  t3 = -1.000000000000001 * coords[0] + coords[4];
  t4 = -0.9999999999999988 * coords[1] + 0.9999999999999988 * coords[3];
  t5 = t1 * t2 + -1.0 * t3 * t4;
  t6 = fabs(t5);
  t8 = 1.0 / t5;
  t9 = t8 * -1.0 * t4;
  t10 = t2 * t8;
  t11 = t1 * t8;
  t12 = t8 * -1.0 * t3;
  for (int32_t ip = 0; ip <= 6; ++ip)
  {
    t18 = 0.0;
    t17 = 0.0;
    t16 = 0.0;
    t15 = 0.0;
    for (int32_t i = 0; i <= 6; ++i)
    {
      t15 = t15 + t14[7 * ip + i] * w_1[i];
      t16 = t16 + t13[7 * ip + i] * w_1[i];
      t17 = t17 + t14[7 * ip + i] * w_0[i];
      t18 = t18 + t13[7 * ip + i] * w_0[i];
    }
    t19 = t12 * t18 + t11 * t17;
    t20 = t12 * t16 + t11 * t15;
    t21 = (-1.0 * t19 * 2.0 * t20 + (t10 * t18 + t9 * t17) * (t10 * t16 + t9 * t15) + t19 * t20) * t7[ip] * t6;
    for (int32_t j = 0; j <= 6; ++j)
      A[j] = A[j] + t0[7 * ip + j] * t21;
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