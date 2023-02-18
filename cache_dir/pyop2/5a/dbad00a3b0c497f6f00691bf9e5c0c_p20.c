#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void expression_kernel(double *__restrict__ C0, double const *__restrict__ C1, double const *__restrict__ C2);
static void expression_kernel(double *__restrict__ C0, double const *__restrict__ C1, double const *__restrict__ C2)
{
  C0[0] = 0.3333333333333333 * C0[0] + 0.6666666666666666 * (C2[0] + C1[0]);
}

void wrap_expression_kernel(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ dat1, double const *__restrict__ dat2)
{
  for (int32_t n = start; n <= -1 + end; ++n)
    expression_kernel(&(dat0[n]), &(dat1[n]), &(dat2[n]));
}