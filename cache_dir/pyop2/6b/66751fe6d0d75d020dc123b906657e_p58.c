#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void form1_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_4);
static void form1_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_4)
{
  double const t0[7 * 7] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
  double t1;
  double t2;
  double const t3[7] = { 0.025, 0.025, 0.025, 0.06666666666666667, 0.06666666666666667, 0.06666666666666667, 0.225 };
  double t4;
  double t5[7];

  t1 = 1.0 / w_4[0];
  t2 = fabs((-0.9999999999999988 * coords[0] + 0.9999999999999988 * coords[2]) * (-1.000000000000001 * coords[1] + coords[5]) + -1.0 * (-1.000000000000001 * coords[0] + coords[4]) * (-0.9999999999999988 * coords[1] + 0.9999999999999988 * coords[3]));
  for (int32_t k0 = 0; k0 <= 6; ++k0)
    t5[k0] = 0.0;
  for (int32_t ip = 0; ip <= 6; ++ip)
  {
    t4 = t1 * t3[ip] * t2;
    for (int32_t k0_0 = 0; k0_0 <= 6; ++k0_0)
      t5[k0_0] = t5[k0_0] + t4 * t0[7 * ip + k0_0] * t0[7 * ip + k0_0];
  }
  for (int32_t k0_1 = 0; k0_1 <= 6; ++k0_1)
  {
    A[2 * k0_1] = A[2 * k0_1] + t5[k0_1];
    A[1 + 2 * k0_1] = A[1 + 2 * k0_1] + t5[k0_1];
  }
}

void wrap_form1_cell_integral_otherwise(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ dat1, double const *__restrict__ glob0, int32_t const *__restrict__ map0, int32_t const *__restrict__ map1)
{
  double t0[7 * 2];
  double t1[3 * 2];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    {
      int32_t const i13 = 0;

      for (int32_t i14 = 0; i14 <= 2; ++i14)
        for (int32_t i15 = 0; i15 <= 1; ++i15)
          t1[2 * i14 + i15] = dat1[2 * map1[3 * n + i14] + i15];
    }
    {
      int32_t const i10 = 0;

      for (int32_t i11 = 0; i11 <= 6; ++i11)
        for (int32_t i12 = 0; i12 <= 1; ++i12)
          t0[2 * i11 + i12] = 0.0;
    }
    form1_cell_integral_otherwise(&(t0[0]), &(t1[0]), &(glob0[0]));
    for (int32_t i8 = 0; i8 <= 1; ++i8)
    {
      int32_t const i9 = 0;

      for (int32_t i7 = 0; i7 <= 6; ++i7)
        dat0[2 * map0[7 * n + i7] + i8] = dat0[2 * map0[7 * n + i7] + i8] + t0[2 * i7 + i8];
    }
  }
}