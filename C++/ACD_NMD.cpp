#include "mex.h"
#include <stdlib.h>
#include<math.h>
#include <algorithm>
#include "GSL_utils.h"
#include <chrono>

double findmin(double* a0, double* b0, double* c0, int m, double* a, double* b, double* c, double* sga, double* aa, double* bb, double* cc, double* ab, double* ac, double* bc, double* breakpts);
void updateH(double* X, double* W, double* H, int m, int n, int r, double* WH, double* aF, double* bF, double* cF, double* a, double* b, double* c, double* sga, double* aa, double* bb, double* cc, double* ab, double* ac, double* bc, double* breakpts);
void updateW(double* X, double* W, double* H, int m, int n, int r, double* WH, double* aF, double* bF, double* cF, double* a, double* b, double* c, double* sga, double* aa, double* bb, double* cc, double* ab, double* ac, double* bc, double* breakpts);
void extrapolation(gsl_matrix* A, gsl_matrix* B, gsl_matrix* C, double a);



//**************************************************************************************//
//																						//
//			Main Function called in Matlab to Compute an approximate solution			//
//			of the following non linear matrix decomposition problem (NMD)				//
//					This algorithm is an alternate scheme for NMD						//
//						using extrapolation for the acceleration						//
//																						//
//				min_{X,W,H} ||X - max(0,WH)||_F  s.t. rank(WH)<=r						//
//																						//
// Using block coordinate descent														//
//																						//
// Update of each element of H															//
// Extrapolation of H																	//
// Update of each element of W															//
// Extrapolation of W																	//
//																						//
// -------------------------------------- Inputs --------------------------------------	//
//					-  X 		: matrix (m x n)										//
//					-  W		: matrix (m x r)										//
//					-  H		: matrix (r x n)										//
//					- maxiter	: the number of iterations for the algorithm			//
//					-  WH 		: matrix (m x n) product of W and H						//
//					- normX		: Frobenius norm of the matrix X						//
//					- 4 Coefficients for extrapolation :								//
//						- beta	(between 0 and 1)										//
//						- eta		( > 1 )												//
//						- gamma		( > 1 )												//
// 						- gamma0	( > 1 )												//
//					- e			: the first error after the scaling						//
//					- tol		: tolerance on the error								//
//																						//
// -------------------------------------- Outputs -------------------------------------	//
//					- Hend		: The new matrix H										//
//					- Wend		: The new matrix W										//
//					- err		: The vector of errors  ||X - max(0,WH)||_F / ||X||_F	//
//																						//
//**************************************************************************************//
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	//Get the input data
	double* X = mxGetPr(prhs[0]);
	double* W = mxGetPr(prhs[1]);
	double* H = mxGetPr(prhs[2]);
	double* maxiter = mxGetPr(prhs[3]);
	double* WH = mxGetPr(prhs[4]);
	double* normX = mxGetPr(prhs[5]);
	double* beta = mxGetPr(prhs[6]);
	double* eta = mxGetPr(prhs[7]);
	double* gamma = mxGetPr(prhs[8]);
	double* gamma0 = mxGetPr(prhs[9]);
	double* e = mxGetPr(prhs[10]);
	double* tol = mxGetPr(prhs[11]);

	//Parameters of extrapolation
	double beta0 = 1.0;
	double beta_prec = *beta;

	//Get the size of the matrix
	int m = mxGetM(prhs[1]);
	int r = mxGetN(prhs[1]);
	int n = mxGetN(prhs[2]);

	//Create the ouput data
	plhs[0] = mxCreateDoubleMatrix(r, n, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(m, r, mxREAL);
	plhs[2] = mxCreateDoubleMatrix((int)*maxiter + 1, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix((int)*maxiter, 1, mxREAL);
	double* Hend = mxGetPr(plhs[0]);
	double* Wend = mxGetPr(plhs[1]);
	double* err = mxGetPr(plhs[2]);
	double* time = mxGetPr(plhs[3]);

	//Initialization of every vectors used in the code
	int size;
	if (n < m) size = m;
	else size = n;
	double* aH = (double*)malloc(m * sizeof(double));
	double* bH = (double*)malloc(m * sizeof(double));
	double* cH = (double*)malloc(m * sizeof(double));
	double* aW = (double*)malloc(n * sizeof(double));
	double* bW = (double*)malloc(n * sizeof(double));
	double* cW = (double*)malloc(n * sizeof(double));
	double* aF = (double*)malloc(size * sizeof(double));
	double* bF = (double*)malloc(size * sizeof(double));
	double* cF = (double*)malloc(size * sizeof(double));
	double* sga = (double*)malloc(size * sizeof(double));
	double* aa = (double*)malloc(size * sizeof(double));
	double* bb = (double*)malloc(size * sizeof(double));
	double* cc = (double*)malloc(size * sizeof(double));
	double* ab = (double*)malloc(size * sizeof(double));
	double* ac = (double*)malloc(size * sizeof(double));
	double* bc = (double*)malloc(size * sizeof(double));
	double* breakpts = (double*)malloc(size * sizeof(double));

	//Conversion of Matrix to GSL format for extrapolation
	gsl_matrix* H_old = mxmat_to_gslmat(prhs[2]);
	gsl_matrix* H_extra = gsl_matrix_alloc(r, n);
	gsl_matrix* H_new = gsl_matrix_alloc(r, n);
	gsl_matrix* W_old = mxmat_to_gslmat(prhs[1]);
	gsl_matrix* W_new = gsl_matrix_alloc(m, r);
	gsl_matrix* W_extra = gsl_matrix_alloc(m, r);
	gsl_matrix* WH_old = mxmat_to_gslmat(prhs[4]);
	gsl_matrix* WH_extra = gsl_matrix_alloc(m, n);

	//Update H, W and calulate the error and the time
	err[0] = *e;

	clock_t start_time = clock();
	clock_t end_time;

	for (int i = 0; i < (int)*maxiter; i++) {

		//Update H
		updateH(X, W, H, m, n, r, WH, aF, bF, cF, aH, bH, cH, sga, aa, bb, cc, ab, ac, bc, breakpts);

		//Extrapolation of H
		//H_new = double_to_gslmat(H, r, n);
		//extrapolation(H_old, H_new, H_extra, *beta);
		//gslmat_to_double(H_extra, H);

		//Matrix Product WH with the new matrix H extrapolated
		//gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, W_old, H_extra, 0.0, WH_extra);
		//gslmat_to_double(WH_extra, WH);

		//Update W
		updateW(X, W, H, m, n, r, WH, aF, bF, cF, aW, bW, cW, sga, aa, bb, cc, ab, ac, bc, breakpts);

		//Extrapolation of W
		//W_new = double_to_gslmat(W, m, r);
		//extrapolation(W_old, W_new, W_extra, *beta);
		//gslmat_to_double(W_extra, W);

		//Matrix Product WH with the new matrix W extrapolated
		//gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, W_extra, H_extra, 0.0, WH_extra);
		//gslmat_to_double(WH_extra, WH);

		//Compute error
		*e = 0.0;
		for (int j = 0; j < m; j++) {
			for (int k = 0; k < n; k++) {
				if (WH[j + k * m] > 0)
					*e += pow(X[j + k * m] - WH[j + k * m], 2);
				else
					*e += pow(X[j + k * m], 2);
			}
		}
		*e = sqrt(*e);
		err[i + 1] = *e / (*normX);

		/*if (err[i + 1] < err[i]) {
			//Change parameters of extrapolation
			*beta = (beta0 < *gamma* (*beta) ? beta0 : *gamma * (*beta));
			beta0 = (1 < *gamma0 * beta0 ? 1 : *gamma0 * beta0);
			beta_prec = *beta;

			//Change the old H, W and WH to the matrix extrapolated
			gsl_matrix_memcpy(H_old, H_extra);
			gsl_matrix_memcpy(W_old, W_extra);
			gsl_matrix_memcpy(WH_old, WH_extra);
		}
		else {
			//Change parameters of extrapolation
			*beta = *beta / (*eta);
			beta0 = beta_prec;

			//Restart H, W and WH to the old values
			gslmat_to_double(H_old, H);
			gslmat_to_double(W_old, W);
			gslmat_to_double(WH_old, WH);
		}*/

		//Compute the time
		end_time = clock();
		time[i] = static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC;

		//Break
		if (err[i + 1] < *tol)
			break;

	}

	//Copy W for output
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < r; j++) {
			Wend[i + j * m] = W[i + j * m];
		}
	}
	//Copy H for output
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < n; j++) {
			Hend[i + j * r] = H[i + j * r];
		}
	}

	//Clean memory
	free(aH);
	free(bH);
	free(cH);
	free(aW);
	free(bW);
	free(cW);
	free(aF);
	free(bF);
	free(cF);
	free(sga);
	free(aa);
	free(bb);
	free(cc);
	free(ab);
	free(ac);
	free(bc);
	free(breakpts);
	gsl_matrix_free(H_old);
	gsl_matrix_free(H_new);
	gsl_matrix_free(H_extra);
	gsl_matrix_free(W_old);
	gsl_matrix_free(W_new);
	gsl_matrix_free(W_extra);
	gsl_matrix_free(WH_old);
	gsl_matrix_free(WH_extra);
}

/****************************************************************************************/
//					 Function to extrapolate a matrix after and update					//
//							A_extra = A_new + beta*(A_new - A_old);						//
//			Inputs : - 3 matrices A_old, A_new and A_extra  (with the same size)		//
//					 - double beta : coefficient of extrapolation						//
/****************************************************************************************/
//void extrapolation(gsl_matrix* A_old, gsl_matrix* A_new, gsl_matrix* A_extra, double beta) {
//	gsl_matrix_memcpy(A_extra, A_new);
//	gsl_matrix_sub(A_extra, A_old);
//	gsl_matrix_scale(A_extra, beta);
//	gsl_matrix_add(A_extra, A_new);
//}

/****************************************************************************************/
//							 Function to update the W matrix							//
//			Inputs : - 3 matrices X, W and H  (pointers double*)						//
//					 - int m, n, r : the dimensions of the matrices						//
//					 - 1 matrix WH, the product W*H that is conserved  					//
//					 - double* err : the error ||X - max(0,WH)||_F / ||X||_F 			//
/****************************************************************************************/
void updateW(double* X, double* W, double* H, int m, int n, int r, double* WH, double* aF, double* bF, double* cF, double* a, double* b, double* c, double* sga, double* aa, double* bb, double* cc, double* ab, double* ac, double* bc, double* breakpts) {

	// Initialization of temporary pointers
	double Wijold;

	// Loop to iterate on each element of H to update them
	// Originally, we called the updateH function with X', H' and W' as parameters
	// But we created this function to avoid computing the transposes
	for (int j = 0; j < m; j++) {

		// c = X(j,:)
		for (int k = 0; k < n; k++) {
			c[k] = X[j + k * m];
		}

		for (int i = 0; i < r; i++) {

			// Saving the old value of W(j,i)
			Wijold = W[j + i * m];

			// a = H(i,:) and b = WH(j,:) - H(i,:)*W(j,i)
			for (int k = 0; k < n; k++) {
				a[k] = H[i + k * r];
				b[k] = WH[j + k * m] - H[i + k * r] * Wijold;
			}

			// call to the function to find the best new value for W(j,i)
			W[j + i * m] = findmin(a, b, c, n, aF, bF, cF, sga, aa, bb, cc, ab, ac, bc, breakpts);

			// WH(j,:) = WH(j,:) + H(i,:)*(W(j,i) - Wijold)
			for (int k = 0; k < n; k++) {
				WH[j + k * m] += H[i + k * r] * (W[j + i * m] - Wijold);
			}
		}
	}
}


/****************************************************************************************/
//							 Function to update the W matrix							//
//			Inputs : - 3 matrices X, W and H  (pointers double*)						//
//					 - int m, n, r : the dimensions of the matrices						//
//					 - 1 matrix WH, the product W*H that is conserved 					//
/****************************************************************************************/
void updateH(double* X, double* W, double* H, int m, int n, int r, double* WH, double* aF, double* bF, double* cF, double* a, double* b, double* c, double* sga, double* aa, double* bb, double* cc, double* ab, double* ac, double* bc, double* breakpts) {

	// Initialization of the old value
	double Hijold;

	// Loop to iterate on each element of H to update them
	for (int i = 0; i < n; i++) {

		// c = X(:,i)
		for (int k = 0; k < m; k++) {
			c[k] = X[k + i * m];
		}

		for (int j = 0; j < r; j++) {
			// Saving the old value of H(i,j)
			Hijold = H[j + i * r];

			// a = W(:,j) and b = WH(:,i) - W(:,j)*H(i,j)
			for (int k = 0; k < m; k++) {
				a[k] = W[k + j * m];
				b[k] = WH[k + i * m] - W[k + j * m] * Hijold;
			}
			// call to the function to find the best new value for H(i,j)
			H[j + i * r] = findmin(a, b, c, m, aF, bF, cF, sga, aa, bb, cc, ab, ac, bc, breakpts);

			// WH(:,i) = WH(:,i) + W(:,j)*(H(i,j) - Hijold)
			for (int k = 0; k < m; k++) {
				WH[k + i * m] += W[k + j * m] * (H[j + i * r] - Hijold);
			}
		}

	}




}


/****************************************************************************************/
//							 Function to optimize the objective							//
//			Inputs : - 3 vectors a,b,c  (pointers double*)								//
//					- int m : the dimension of the vectors								//
//			Output : xopt = argmin_x ||c-max(0,b+ax)||_2^2 (type double)				//
/****************************************************************************************/
double findmin(double* a0, double* b0, double* c0, int m, double* a, double* b, double* c, double* sga, double* aa, double* bb, double* cc, double* ab, double* ac, double* bc, double* breakpts) {

	/***** Delete the elements for the index where a == 0 *****/

	// Get the values of a,b,c for the index where a != 0
	int nnz = 0;
	for (int i = 0; i < m; i++)
	{
		if (abs(a0[i]) > pow(10, -16))
		{
			a0[nnz] = a0[i];
			b0[nnz] = b0[i];
			c0[nnz] = c0[i];
			nnz = nnz + 1;
		}
	}

	/***** Breakpoints computing ******/

	double* breakpts_nnz = (double*)malloc(nnz * sizeof(double));
	for (int i = 0; i < nnz; i++) breakpts_nnz[i] = -b0[i] / a0[i];

	// Breakpoints sorting to iterate in the right order + a, b, c sorting in the same order
	std::size_t* indices = (std::size_t*)malloc(nnz * sizeof(std::size_t));

	//Create a vector of index
	for (int i = 0; i < nnz; i++) {
		indices[i] = i;
	}

	//Get the index when sort the breakpoints
	gsl_sort_index(indices, breakpts_nnz, 1, nnz);

	//put in order the elements of vectors a,b,c and breakpts
	for (int i = 0; i < nnz; i++) {
		breakpts[i] = breakpts_nnz[(int)indices[i]];
		a[i] = a0[(int)indices[i]];
		b[i] = b0[(int)indices[i]];
		c[i] = c0[(int)indices[i]];
	}

	//clean memory
	free(indices);
	free(breakpts_nnz);

	/***** Intermediate calculus for coefficient computing *****/
	for (int i = 0; i < nnz; i++) {
		sga[i] = (a[i] < 0 ? -1 : 1);
		aa[i] = a[i] * a[i];
		bb[i] = b[i] * b[i];
		cc[i] = c[i] * c[i];
		ab[i] = a[i] * b[i];
		ac[i] = a[i] * c[i];
		bc[i] = b[i] * c[i];
	}

	/***** Quadratic function far left *****/
	double ti = 0.0; double tl = 0.0; double tq = 0.0;

	//Calculating coefficients for the quadratic function where a is negative
	//Because we dont need the constant
	for (int i = 0; i < nnz; i++)
	{
		ti += cc[i];
		if (a[i] < pow(-10, -16))
		{
			ti = ti - 2 * bc[i] + bb[i];
			tl = tl + 2 * (ab[i] - ac[i]);
			tq = tq + aa[i];
		}
	}

	//The minimum of tq*x^2+tl*x+ti is at xmin=-tl/(2*tq)
	double xmin = -tl / (2 * tq);

	//Which is the best: the min or the first breakpoint ?
	double xopt = (xmin < breakpts[0] ? xmin : breakpts[0]);

	//The value of f at xopt
	double yopt = ti + tl * (xopt)+tq * pow(xopt, 2);

	//For each interval, we check if the min (or the previous breakpoint) is better
	double xt;
	for (int i = 0; i < nnz; i++) {

		//Update the coefficients of the quadratic function
		ti += sga[i] * (bb[i] - 2 * bc[i]);
		tl += sga[i] * (2 * ab[i] - 2 * ac[i]);
		tq += sga[i] * (aa[i]);

		//the minimum of tq*x^2+tl*x+ti is at xmin=-tl/(2*tq)
		xmin = -(tl) / (2 * tq);

		//Check if xmin is between two breakpoints
		if (((i + 1 < nnz) && ((xmin < breakpts[i + 1]) && (xmin > breakpts[i]))) || ((i + 1 >= nnz) && (xmin > breakpts[nnz - 1])))
			xt = xmin;
		else
			xt = breakpts[i];

		//Check if the f(xt) is lower than yopt
		if (yopt > ti + tl * xt + tq * pow(xt, 2)) {
			xopt = xt;
			yopt = ti + tl * (xopt)+tq * pow(xopt, 2);
		}

	}

	return xopt;

}

