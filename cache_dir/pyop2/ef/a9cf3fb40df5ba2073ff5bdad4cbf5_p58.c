#include <evaluate.h>
#include <math.h>
struct ReferenceCoords {
    double X[2];
};

static inline void to_reference_coords_kernel(void *result_, double *x0, int *return_value, double *C)
{
    struct ReferenceCoords *result = (struct ReferenceCoords *) result_;

    /*
     * Mapping coordinates from physical to reference space
     */

    double *X = result->X;
    X[0] = 0.3;
X[1] = 0.3;

    int converged = 0;
    for (int it = 0; !converged && it < 1; it++) {
        double dX[2] = { 0.0 };
double  t0  = ((((((-0.999999999999999) * (C[0]))) + (((0.999999999999999) * (C[2]))))) + (((0) * (C[4]))));
double  t1  = ((((((-1) * (C[1]))) + (((0) * (C[3]))))) + (C[5]));
double  t2  = ((((((-1) * (C[0]))) + (((0) * (C[2]))))) + (C[4]));
double  t3  = ((((((-0.999999999999999) * (C[1]))) + (((0.999999999999999) * (C[3]))))) + (((0) * (C[5]))));
double  t4  = ((((t0) * (t1))) + (((-1) * (((t2) * (t3))))));
double  t5  = ((((-1.73205080756888) + (((1.73205080756888) * (X[1]))))) + (((3.46410161513775) * (X[0]))));
double  t6  = ((-1) + (((3) * (X[1]))));
double  t7  = ((-0.166666666666667) * (t6));
double  t8  = ((((0.333333333333333) + (((-0.288675134594813) * (t5))))) + (t7));
double  t9  = ((((0.333333333333333) + (((0.288675134594813) * (t5))))) + (t7));
double  t10  = ((((0.333333333333333) + (((0) * (t5))))) + (((0.333333333333333) * (t6))));
double  t11  = ((((((((C[0]) * (t8))) + (((C[2]) * (t9))))) + (((C[4]) * (t10))))) + (((-1) * (x0[0]))));
double  t12  = ((((((((C[1]) * (t8))) + (((C[3]) * (t9))))) + (((C[5]) * (t10))))) + (((-1) * (x0[1]))));
dX[0] += (((((t1) / (t4))) * (t11))) + (((((((-1) * (t2))) / (t4))) * (t12)));
dX[1] += (((((((-1) * (t3))) / (t4))) * (t11))) + (((((t0) / (t4))) * (t12)));


        if (PetscRealPart(dX[0])*PetscRealPart(dX[0]) + PetscRealPart(dX[1])*PetscRealPart(dX[1]) < 1e-12 * 1e-12) {
            converged = 1;
        }

	X[0] -= dX[0];
	X[1] -= dX[1];
    }

    // Are we inside the reference element?
    *return_value = (PetscRealPart(X[0]) + 0.01 >= 0) && (PetscRealPart(X[1]) + 0.01 >= 0) && (PetscRealPart(X[0]) + PetscRealPart(X[1]) - 0.01 <= 1);
}

static inline void wrap_to_reference_coords(
    void* const result_, double* const x, int* const return_value, int32_t const start, int32_t const end,
    double const *__restrict__ coords, int32_t const *__restrict__ coords_map);

int to_reference_coords(void *result_, struct Function *f, int cell, double *x)
{
    int return_value = 0;
    wrap_to_reference_coords(result_, x, &return_value, cell, cell+1, f->coords, f->coords_map);
    return return_value;
}

int to_reference_coords_xtr(void *result_, struct Function *f, int cell, int layer, double *x)
{
    int return_value = 0;
    //int layers[2] = {0, layer+2};  // +2 because the layer loop goes to layers[1]-1, which is nlayers-1
    //wrap_to_reference_coords(result_, x, &return_value, cell, cell+1, layers, f->coords, f->coords_map);
    return return_value;
}


#include <complex.h>
#include <math.h>
#include <petsc.h>

#include <stdint.h>

void wrap_to_reference_coords(void* const farg0, double* const farg1, int* const farg2, int32_t const start, int32_t const end, double const *__restrict__ dat0, int32_t const *__restrict__ map0)
{
  double t0[3 * 2];

  {
    int32_t const n = start;

    {
      int32_t const i3 = 0;

      for (int32_t i4 = 0; i4 <= 2; ++i4)
        for (int32_t i5 = 0; i5 <= 1; ++i5)
          t0[2 * i4 + i5] = dat0[2 * map0[3 * start + i4] + i5];
    }
    to_reference_coords_kernel(farg0, farg1, farg2, &(t0[0]));
  }
}
#include <stdio.h>
#include <stdlib.h>
#include <spatialindex/capi/sidx_api.h>

#include <evaluate.h>

int locate_cell(struct Function *f,
        double *x,
        int dim,
        inside_predicate try_candidate,
        inside_predicate_xtr try_candidate_xtr,
        void *data_)
{
    RTError err;
    int cell = -1;

    if (f->sidx) {
        int64_t *ids = NULL;
        uint64_t nids = 0;
        err = Index_Intersects_id(f->sidx, x, x, dim, &ids, &nids);
        if (err != RT_None) {
            fputs("ERROR: Index_Intersects_id failed in libspatialindex!", stderr);
            return -1;
        }
        if (f->extruded == 0) {
            for (uint64_t i = 0; i < nids; i++) {
                if ((*try_candidate)(data_, f, ids[i], x)) {
                    cell = ids[i];
                    break;
                }
            }
        }
        else {
            for (uint64_t i = 0; i < nids; i++) {
                int nlayers = f->n_layers;
                int c = ids[i] / nlayers;
                int l = ids[i] % nlayers;
                if ((*try_candidate_xtr)(data_, f, c, l, x)) {
                    cell = ids[i];
                    break;
                }
            }
        }
        free(ids);
    } else {
        if (f->extruded == 0) {
            for (int c = 0; c < f->n_cols; c++) {
                if ((*try_candidate)(data_, f, c, x)) {
                    cell = c;
                    break;
                }
            }
        }
        else {
            for (int c = 0; c < f->n_cols; c++) {
                for (int l = 0; l < f->n_layers; l++)
                    if ((*try_candidate_xtr)(data_, f, c, l, x)) {
                        cell = l;
                        break;
                    }
                if (cell != -1) {
                    cell = c * f->n_layers + cell;
                    break;
                }
            }
        }
    }
    return cell;
}

    int locator(struct Function *f, double *x, double *X)
    {
        struct ReferenceCoords reference_coords;
        int cell = locate_cell(f, x, 2, &to_reference_coords, &to_reference_coords_xtr, &reference_coords);
        for(int i=0; i<2; i++) {
            X[i] = reference_coords.X[i];
        }
        return cell;
    }
    