#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void form0_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0);
static void form0_cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0)
{
  double t0;
  double const t1[1] = { 1.0 };

  t0 = (w_0[0] < 1.5 ? 1.0 : 0.0) * 0.5 * fabs((-0.9999999999999988 * coords[0] + 0.9999999999999988 * coords[2]) * (-1.0 * coords[1] + coords[5]) + -1.0 * (-1.0 * coords[0] + coords[4]) * (-0.9999999999999988 * coords[1] + 0.9999999999999988 * coords[3]));
  {
    int32_t const j = 0;

    A[0] = A[0] + t1[0] * t0;
  }
}

void wrap_form0_cell_integral_otherwise(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ dat1, double const *__restrict__ dat2, int32_t const *__restrict__ map0, int32_t const *__restrict__ map1)
{
  double t0[1];
  double t1[3 * 2];
  double t2[1];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    {
      int32_t const i18 = 0;

      {
        int32_t const i19 = 0;

        {
          int32_t const i20 = 0;

          t2[0] = dat2[map0[n]];
        }
      }
    }
    {
      int32_t const i15 = 0;

      for (int32_t i16 = 0; i16 <= 2; ++i16)
        for (int32_t i17 = 0; i17 <= 1; ++i17)
          t1[2 * i16 + i17] = dat1[2 * map1[3 * n + i16] + i17];
    }
    {
      int32_t const i12 = 0;

      {
        int32_t const i13 = 0;

        {
          int32_t const i14 = 0;

          t0[0] = 0.0;
        }
      }
    }
    form0_cell_integral_otherwise(&(t0[0]), &(t1[0]), &(t2[0]));
    {
      int32_t const i10 = 0;

      {
        int32_t const i9 = 0;

        {
          int32_t const i11 = 0;

          dat0[map0[n]] = dat0[map0[n]] + t0[0];
        }
      }
    }
  }
}