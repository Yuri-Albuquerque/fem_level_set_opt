#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void form1_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords);
static void form1_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords)
{
  double const t0[7 * 7] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
  double t1;
  double const t2[7] = { 0.025, 0.025, 0.025, 0.06666666666666667, 0.06666666666666667, 0.06666666666666667, 0.225 };
  double t3;
  double t4[7];

  t1 = fabs((-0.9999999999999988 * coords[0] + 0.9999999999999988 * coords[2]) * (-1.000000000000001 * coords[1] + coords[5]) + -1.0 * (-1.000000000000001 * coords[0] + coords[4]) * (-0.9999999999999988 * coords[1] + 0.9999999999999988 * coords[3]));
  for (int32_t k0 = 0; k0 <= 6; ++k0)
    t4[k0] = 0.0;
  for (int32_t ip = 0; ip <= 6; ++ip)
  {
    t3 = 10000.0 * t2[ip] * t1;
    for (int32_t k0_0 = 0; k0_0 <= 6; ++k0_0)
      t4[k0_0] = t4[k0_0] + t3 * t0[7 * ip + k0_0] * t0[7 * ip + k0_0];
  }
  for (int32_t k0_1 = 0; k0_1 <= 6; ++k0_1)
  {
    A[2 * k0_1] = A[2 * k0_1] + t4[k0_1];
    A[1 + 2 * k0_1] = A[1 + 2 * k0_1] + t4[k0_1];
  }
}

void wrap_form1_cell_integral_otherwise(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ dat1, int32_t const *__restrict__ map0, int32_t const *__restrict__ map1)
{
  double t0[7 * 2];
  double t1[3 * 2];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    {
      int32_t const i9 = 0;

      for (int32_t i10 = 0; i10 <= 6; ++i10)
        for (int32_t i11 = 0; i11 <= 1; ++i11)
          t0[2 * i10 + i11] = 0.0;
    }
    {
      int32_t const i12 = 0;

      for (int32_t i13 = 0; i13 <= 2; ++i13)
        for (int32_t i14 = 0; i14 <= 1; ++i14)
          t1[2 * i13 + i14] = dat1[2 * map1[3 * n + i13] + i14];
    }
    form1_cell_integral_otherwise(&(t0[0]), &(t1[0]));
    for (int32_t i7 = 0; i7 <= 1; ++i7)
    {
      int32_t const i8 = 0;

      for (int32_t i6 = 0; i6 <= 6; ++i6)
        dat0[2 * map0[7 * n + i6] + i7] = dat0[2 * map0[7 * n + i6] + i7] + t0[2 * i6 + i7];
    }
  }
}