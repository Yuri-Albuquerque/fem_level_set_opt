#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void expression_kernel(double *__restrict__ C0, double const *__restrict__ C1, double const *__restrict__ C2, double const *__restrict__ C3);
static void expression_kernel(double *__restrict__ C0, double const *__restrict__ C1, double const *__restrict__ C2, double const *__restrict__ C3)
{
  C0[0] = 0.75 * C3[0] + 0.25 * (C2[0] + C1[0]);
}

void wrap_expression_kernel(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ dat1, double const *__restrict__ dat2, double const *__restrict__ dat3)
{
  for (int32_t n = start; n <= -1 + end; ++n)
    expression_kernel(&(dat0[n]), &(dat1[n]), &(dat2[n]), &(dat3[n]));
}