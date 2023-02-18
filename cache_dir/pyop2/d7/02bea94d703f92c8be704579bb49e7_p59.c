#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void expression_kernel(double *__restrict__ C0);
static void expression_kernel(double *__restrict__ C0)
{
  C0[0] = 4.4;
}

void wrap_expression_kernel(int32_t const start, int32_t const end, double *__restrict__ dat0)
{
  for (int32_t n = start; n <= -1 + end; ++n)
    expression_kernel(&(dat0[n]));
}