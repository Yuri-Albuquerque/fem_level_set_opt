#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void form0_exterior_facet_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2, uint32_t const *__restrict__ facet);
static void form0_exterior_facet_integral_otherwise(double *__restrict__ A, double const *__restrict__ coords, double const *__restrict__ w_0, double const *__restrict__ w_1, double const *__restrict__ w_2, uint32_t const *__restrict__ facet)
{
  double const t0[3 * 2] = { -1.0, 1.0, 0.0, 1.0, 1.0, 0.0 };
  double t1;
  double t2;
  double t3;
  double const t4[4] = { 0.17392742256872684, 0.3260725774312731, 0.3260725774312731, 0.17392742256872684 };
  double const t5[3 * 4 * 7] = { 0.0, 0.8013460293699309, -0.05979028222412171, 0.25844425285419076, 0.0, 0.0, 0.0, 0.0, 0.227784076790952, -0.11219696679390412, 0.884412890002952, 0.0, 0.0, 0.0, 0.0, -0.11219696679390408, 0.227784076790952, 0.8844128900029522, 0.0, 0.0, 0.0, 0.0, -0.059790282224121666, 0.8013460293699308, 0.2584442528541909, 0.0, 0.0, 0.0, 0.8013460293699309, 0.0, -0.05979028222412171, 0.0, 0.2584442528541907, 0.0, 0.0, 0.227784076790952, 0.0, -0.11219696679390412, 0.0, 0.884412890002952, 0.0, 0.0, -0.11219696679390412, 0.0, 0.2277840767909521, 0.0, 0.8844128900029522, 0.0, 0.0, -0.059790282224121666, 0.0, 0.8013460293699308, 0.0, 0.25844425285419087, 0.0, 0.0, 0.8013460293699308, -0.05979028222412166, 0.0, 0.0, 0.0, 0.2584442528541908, 0.0, 0.22778407679095206, -0.11219696679390416, 0.0, 0.0, 0.0, 0.8844128900029519, 0.0, -0.11219696679390415, 0.22778407679095206, 0.0, 0.0, 0.0, 0.8844128900029521, 0.0, -0.059790282224121714, 0.8013460293699307, 0.0, 0.0, 0.0, 0.25844425285419104, 0.0 };
  double t6;
  double t7;
  double t8;
  double t9;

  t1 = t0[2 * facet[0]] * (-0.9999999999999988 * coords[0] + 0.9999999999999988 * coords[2]) + t0[2 * facet[0] + 1] * (-1.0 * coords[0] + coords[4]);
  t2 = t0[2 * facet[0]] * (-0.9999999999999988 * coords[1] + 0.9999999999999988 * coords[3]) + t0[2 * facet[0] + 1] * (-1.0 * coords[1] + coords[5]);
  t3 = sqrt(t1 * t1 + t2 * t2);
  for (int32_t ip = 0; ip <= 3; ++ip)
  {
    t8 = 0.0;
    t7 = 0.0;
    t6 = 0.0;
    for (int32_t i = 0; i <= 6; ++i)
    {
      t6 = t6 + t5[28 * facet[0] + 7 * ip + i] * w_0[i];
      t7 = t7 + t5[28 * facet[0] + 7 * ip + i] * w_1[i];
      t8 = t8 + t5[28 * facet[0] + 7 * ip + i] * w_2[i];
    }
    t9 = -1.0 * t6 * (-1.0 * t8 + t7) * 10000.0 * t4[ip] * t3;
    for (int32_t j = 0; j <= 6; ++j)
      A[j] = A[j] + t5[28 * facet[0] + 7 * ip + j] * t9;
  }
}

void wrap_form0_exterior_facet_integral_otherwise(int32_t const start, int32_t const end, double *__restrict__ dat1, double const *__restrict__ dat2, double const *__restrict__ dat3, double const *__restrict__ dat4, double const *__restrict__ dat5, uint32_t const *__restrict__ dat0, int32_t const *__restrict__ map0, int32_t const *__restrict__ map1)
{
  double t0[7];
  double t1[3 * 2];
  double t2[7];
  double t3[7];
  double t4[7];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    for (int32_t i26 = 0; i26 <= 6; ++i26)
    {
      {
        int32_t const i30 = 0;

        {
          int32_t const i31 = 0;

          t4[i26] = dat5[map0[7 * n + i26]];
        }
      }
      {
        int32_t const i28 = 0;

        {
          int32_t const i29 = 0;

          t3[i26] = dat4[map0[7 * n + i26]];
        }
      }
      {
        int32_t const i25 = 0;

        {
          int32_t const i27 = 0;

          t2[i26] = dat3[map0[7 * n + i26]];
        }
      }
    }
    {
      int32_t const i22 = 0;

      for (int32_t i23 = 0; i23 <= 2; ++i23)
        for (int32_t i24 = 0; i24 <= 1; ++i24)
          t1[2 * i23 + i24] = dat2[2 * map1[3 * n + i23] + i24];
    }
    {
      int32_t const i19 = 0;

      for (int32_t i20 = 0; i20 <= 6; ++i20)
      {
        int32_t const i21 = 0;

        t0[i20] = 0.0;
      }
    }
    form0_exterior_facet_integral_otherwise(&(t0[0]), &(t1[0]), &(t2[0]), &(t3[0]), &(t4[0]), &(dat0[n]));
    {
      int32_t const i17 = 0;

      {
        int32_t const i18 = 0;

        for (int32_t i16 = 0; i16 <= 6; ++i16)
          dat1[map0[7 * n + i16]] = dat1[map0[7 * n + i16]] + t0[i16];
      }
    }
  }
}