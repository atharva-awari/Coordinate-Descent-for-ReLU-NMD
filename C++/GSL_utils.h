#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

//**************************************************************************************//
//																						//
//			This file contains all the functions that we used to convert a matrix		//
//                               or a vector in gsl format                              //
//																						//
//**************************************************************************************//
// convert mx array to GSL matrix
gsl_matrix *mxmat_to_gslmat (const mxArray * mxmat)
{
  const mwSize *dims = mxGetDimensions (mxmat);
  int nb_row = (int) dims[0];
  int nb_col = (int) dims[1];
  double *rawmat = mxGetPr (mxmat);

  gsl_matrix *m = gsl_matrix_alloc (nb_row, nb_col);
  for (int i = 0; i < nb_row; ++i)
    for (int j = 0; j < nb_col; ++j)
      gsl_matrix_set (m, i, j, rawmat[i + j * nb_row]);

  return m;
}

// convert GSL matrix to  mx array
mxArray *gslmat_to_mxmat (const gsl_matrix * gslmat)
{
  int nb_row = gslmat->size1;
  int nb_col = gslmat->size2;

  mxArray *m = mxCreateDoubleMatrix (nb_row, nb_col, mxREAL);
  double *pm = mxGetPr (m);
  for (int i = 0; i < nb_row; ++i)
    for (int j = 0; j < nb_col; ++j)
      pm[i + j * nb_row] = gsl_matrix_get (gslmat, i, j);

  return m;
}

// convert mx array to GSL vector
gsl_vector *mxmat_to_gslvec (const mxArray * mxmat)
{
  const mwSize *dims = mxGetDimensions (mxmat);
  int size = (int) dims[0];
  double *rawvec = mxGetPr (mxmat);

  gsl_vector *m = gsl_vector_alloc (size);
  for (int i = 0; i < size; ++i)
      gsl_vector_set (m, i, rawvec[i]);

  return m;
}

// convert GSL vector to mx array
mxArray *gslvec_to_mxmat (const gsl_vector * gslvec)
{
  int size = gslvec->size;

  mxArray *m = mxCreateDoubleMatrix (size, 1, mxREAL);
  double *pm = mxGetPr (m);
  for (int i = 0; i < size; ++i)
      pm[i] = gsl_vector_get (gslvec, i);

  return m;
}

// convert double vector to GSL vector
gsl_vector *double_to_gslvec(double *array, size_t length) {
    gsl_vector *vector = gsl_vector_alloc(length);

    for (size_t i = 0; i < length; i++) {
        gsl_vector_set(vector, i, array[i]);
    }

    return vector;
}

// convert double vector to GSL matrix
gsl_matrix *double_to_gslmat(double *array, size_t m, size_t n) {
	gsl_matrix *matrix = gsl_matrix_alloc(m, n);

    for (int j=0; j<n; j++) {
        for (int i=0; i<m; i++) {
        	gsl_matrix_set(matrix, i, j, array[i + j*m]);
   		}
    }

    return matrix;
}

// convert GSL vector to double vector
double *gslvec_to_double(gsl_vector *gslvec) {
	int n = gslvec->size;
    double* array = (double*)malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) {
        array[i] = gsl_vector_get(gslvec, i);
    }

    return array;
}

// convert GSL matrix to double vector
void gslmat_to_double(gsl_matrix* gslmat, double* array) {
    int m = gslmat->size1;
    int n = gslmat->size2;

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            array[i + j * m] = gsl_matrix_get(gslmat, i, j);
        }
    }
}
