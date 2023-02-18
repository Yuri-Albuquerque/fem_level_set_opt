#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <petsc.h>
#include <stdint.h>

static void slate_loopy(double *__restrict__ output, double const *__restrict__ T0);
static void slate_wrapper(double *__restrict__ output, double const *__restrict__ coords);
static void subkernel0_cell_to__cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords);
static void slate_loopy(double *__restrict__ output, double const *__restrict__ T0)
{
  {
    int32_t const i = 0;

    {
      int32_t const i_0 = 0;

      output[0] = output[0] + 1.0 / T0[0];
    }
  }
}

static void slate_wrapper(double *__restrict__ output, double const *__restrict__ coords)
{
  double T0[1];

  {
    int32_t const i_id_0 = 0;

    {
      int32_t const i_id = 0;

      T0[0] = 0.0;
    }
  }
  subkernel0_cell_to__cell_integral_otherwise(&(T0[0]), &(coords[0]));
  slate_loopy(&(output[0]), &(T0[0]));
}

static void subkernel0_cell_to__cell_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords)
{
  double const t0[1] = { 1.0 };
  double t1;

  t1 = 0.5 * fabs((-0.9999999999999988 * coords[0] + 0.9999999999999988 * coords[2]) * (-1.0 * coords[1] + coords[5]) + -1.0 * (-1.0 * coords[0] + coords[4]) * (-0.9999999999999988 * coords[1] + 0.9999999999999988 * coords[3]));
  {
    int32_t const j = 0;

    {
      int32_t const k = 0;

      A[0] = A[0] + t1 * t0[0] * t0[0];
    }
  }
}

void wrap_slate_wrapper(int32_t const start, int32_t const end, Mat const mat0, double const *__restrict__ dat0, int32_t const *__restrict__ map0, int32_t const *__restrict__ map1)
{
  double t0[1];
  double t1[3 * 2];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    {
      int32_t const i21 = 0;

      for (int32_t i22 = 0; i22 <= 2; ++i22)
        for (int32_t i23 = 0; i23 <= 1; ++i23)
          t1[2 * i22 + i23] = dat0[2 * map1[3 * n + i22] + i23];
    }
    {
      int32_t const i15 = 0;

      {
        int32_t const i16 = 0;

        {
          int32_t const i17 = 0;

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
        }
      }
    }
    slate_wrapper(&(t0[0]), &(t1[0]));
    MatSetValuesBlockedLocal(mat0, 1, &(map0[n]), 1, &(map0[n]), &(t0[0]), ADD_VALUES);
  }
}