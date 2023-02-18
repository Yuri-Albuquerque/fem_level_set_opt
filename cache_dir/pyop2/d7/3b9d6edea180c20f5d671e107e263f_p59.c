#include <complex.h>
#include <math.h>
#include <petsc.h>
#include <stdint.h>

static void iop_iadd(double *__restrict__ self, double const *__restrict__ other);
static void iop_iadd(double *__restrict__ self, double const *__restrict__ other)
{
  for (int32_t i = 0; i <= 1; ++i)
    self[i] = self[i] + other[i];
}

void wrap_iop_iadd(int32_t const start, int32_t const end, double *__restrict__ dat0, double const *__restrict__ dat1)
{
  for (int32_t n = start; n <= -1 + end; ++n)
    iop_iadd(&(dat0[2 * n]), &(dat1[2 * n]));
}