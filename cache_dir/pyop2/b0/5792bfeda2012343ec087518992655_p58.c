#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void form_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0);
static void form_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0)
{
  double t0;
  double const t1[6] = { 0.054975871827661, 0.054975871827661, 0.054975871827661, 0.1116907948390055, 0.1116907948390055, 0.1116907948390055 };
  double const t2[6 * 7] = { -0.05425305933911809, 0.5381830903967516, -0.05425305933911873, 0.21701223735647515, -0.04865818211316347, 0.21701223735647193, 0.18495673568170157, -0.05425305933911805, -0.05425305933911868, 0.5381830903967515, 0.21701223735647518, 0.21701223735647188, -0.048658182113163556, 0.1849567356817018, 0.5381830903967495, -0.0542530593391185, -0.05425305933911854, -0.04865818211316399, 0.21701223735647399, 0.21701223735647399, 0.1849567356817034, 0.016286982193490795, -0.020235133084975102, 0.01628698219349079, -0.06514792877396303, 0.5374987861648951, -0.0651479287739631, 0.5804582400810245, 0.01628698219349079, 0.01628698219349079, -0.020235133084975133, -0.06514792877396305, -0.06514792877396312, 0.5374987861648951, 0.5804582400810245, -0.02023513308497509, 0.01628698219349081, 0.01628698219349079, 0.5374987861648952, -0.0651479287739631, -0.0651479287739631, 0.5804582400810245 };
  double t3;

  t0 = fabs((-0.9999999999999988 * coords[0] + 0.9999999999999988 * coords[2]) * (-1.0 * coords[1] + coords[5]) + -1.0 * (-1.0 * coords[0] + coords[4]) * (-0.9999999999999988 * coords[1] + 0.9999999999999988 * coords[3]));
  for (int32_t ip = 0; ip <= 5; ++ip)
  {
    t3 = 0.0;
    for (int32_t i = 0; i <= 6; ++i)
      t3 = t3 + t2[7 * ip + i] * w_0[i];
    A[0] = A[0] + t1[ip] * t0 * t3 * t3;
  }
}

void wrap_form_cell_integral_otherwise(int32_t const start, int32_t const end, double *__restrict__ glob0, double const *__restrict__ dat0, double const *__restrict__ dat1, int32_t const *__restrict__ map0, int32_t const *__restrict__ map1)
{
  double t0[1];
  double t1[3 * 2];
  double t2[7];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    {
      int32_t const i9 = 0;

      for (int32_t i10 = 0; i10 <= 2; ++i10)
        for (int32_t i11 = 0; i11 <= 1; ++i11)
          t1[2 * i10 + i11] = dat0[2 * map0[3 * n + i10] + i11];
    }
    {
      int32_t const i8 = 0;

      t0[0] = 0.0;
    }
    {
      int32_t const i12 = 0;

      for (int32_t i13 = 0; i13 <= 6; ++i13)
      {
        int32_t const i14 = 0;

        t2[i13] = dat1[map1[7 * n + i13]];
      }
    }
    form_cell_integral_otherwise(&(t0[0]), &(t1[0]), &(t2[0]));
    {
      int32_t const i7 = 0;

      glob0[0] = glob0[0] + t0[0];
    }
  }
}