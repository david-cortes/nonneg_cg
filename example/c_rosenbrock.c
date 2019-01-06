/*	Example usage of non-negative conjugate gradient optimizer
	with Rosenbrock's (a.k.a. banana) function.
	Implementation taken from here:
	https://gist.github.com/Bismarrck/772d4ae4a8c5ddcd37b8
*/
#include <stdlib.h>
#include <stdio.h>
#include "nonnegcg.h"

void rosen(double x[], int n, double *f, void *data)
{
	double dx = 0.0;
	double sum = 0.0;
	for (int i = 0; i < n - 1; i ++) {
		dx = x[i + 1] - x[i] * x[i];
		sum += 100.0 * dx * dx;
		dx = 1.0 - x[i];
		sum += dx * dx;
	}
	*f = sum;
}

void rosen_der(double x[], int n, double grad[], void *data)
{
	grad[0] = -400.0 * x[0] * (x[1] - x[0] * x[0]) - 2.0 * (1.0 - x[0]);
	grad[n - 1] = 200.0 * (x[n - 1] - x[n - 2] * x[n - 2]);
	
	double d1 = 0.0;
	double d2 = 0.0;
	double d3 = 0.0;
	
	for (int i = 1; i < n - 1; i ++) {
		d1 = 200.0 * (x[i] - x[i - 1] * x[i - 1]);
		d2 = 400.0 * (x[i + 1] - x[i] * x[i]) * x[i];
		d3 = 2.0 * (1.0 - x[i]);
		grad[i] = d1 - d2 - d3;
	}
}

int main()
{
	double x[] = {1.3, 0.7, 0.8, 1.9, 1.2};
	int n = 4;
	printf("Initial values of x: [ ");
	for (size_t i = 0; i < n; i++){printf("%f ", x[i]);} printf("]\n");


	double fun_val;
	size_t niter, nfeval;

	double tol = 1e-6;
	size_t maxnfeval = 1000;
	size_t maxiter = 100;
	double decr_lnsrch = 0.5;
	double lnsrch_const = 0.01;
	size_t max_ls = 20;
	int extra_nonneg_tol = 0;
	int nthreads = 1;
	int verbose = 1;
	double *buffer = (double*) malloc(sizeof(double) * 4 * 4);

	minimize_nonneg_cg(x, n, &fun_val,
		rosen, rosen_der, NULL, NULL,
		tol, maxnfeval, maxiter, &niter, &nfeval,
		decr_lnsrch, lnsrch_const, max_ls,
		extra_nonneg_tol, buffer, nthreads, verbose);

	free(buffer);

	printf("Optimal values of x: [ ");
	for (size_t i = 0; i < n; i++){printf("%f ", x[i]);} printf("]\n");

	return 0;
}
