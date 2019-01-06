// #include <math.h>
// #include <stdlib.h>
// #include <stddef.h>
// #include <limits.h>
// #ifdef _OPENMP
// 	#include <omp.h>
// #endif

// // #ifdef _FOR_PYTON
// // 	#include "findblas.h"
// // #elif defined(_FOR_R)
// // 	#include <R_ext/BLAS.h>
// // 	#include <R_ext/Print.h>
// // 	double cblas_ddot(int n, double *x, int incx, double *y, int incy)
// // 	{
// // 	  return ddot_(&n, x, &incx, y, &incy);
// // 	}
// // 	void cblas_daxpy(int n, double a, double *x, int incx, double *y, int incy)
// // 	{
// // 	  daxpy_(&n, &a, x, &incx, y, &incy);
// // 	}
// // 	void cblas_dscal(int n, double alpha, double *x, int incx)
// // 	{
// // 	  dscal_(&n, &alpha, x, &incx);
// // 	}
// // 	#define printf Rprintf
// // 	#define fprintf(f, message) REprintf(message)
// // #else
// // 	#include "blasfuns.h"
// // #endif

// // #ifndef _FOR_R
// 	#include <stdio.h> 
// // #endif

// /* make a callback print of python and r, pass more params to the callback */

// /* Aliasing for compiler optimizations */
// #ifndef restrict
// 	#ifdef __restrict /* MSVC doesn't have proper support for these optimizations */
// 		#define restrict __restrict
// 	#else
// 		#define restrict 
// 	#endif
// #endif

// /*	OpenMP < 3.0 (e.g. MSVC as of 2019) does not support parallel for's with unsigned iterators,
// 	and does not support declaring the iterator type in the loop itself */
// #ifdef _OPENMP
// 	#if _OPENMP < 20080101 /* OpenMP < 3.0 */
// 		#define size_t_for 
// 	#else
// 		#define size_t_for size_t
// 	#endif
// #else
// 	#define size_t_for size_t
// #endif

// #ifndef isnan
// 	#ifdef _isnan
// 		#define isnan _isnan
// 	#else
// 		#define isnan(x) ( (x) != (x) )
// 	#endif
// #endif
// #ifndef isinf
// 	#ifdef _finite
// 		#define isinf(x) (!_finite(x))
// 	#else
// 		#define isinf(x) ( (x) >= HUGE_VAL || (x) <= -HUGE_VAL )
// 	#endif
// #endif

// #define get_curr_ix_rotation(ix, n) (  ((ix) == 0) ? 0 : (n)  )
// #define incr_ix_rotation(ix) (  ((ix) == 0)? 1 : 0  )
// #define min(a, b) (  ( (a) <= (b) )? a : b  )
// #define square(x) ( (x) * (x) )

// typedef void fun_eval(double x[], int n, double *f, void *data);
// typedef void grad_eval(double x[], int n, double grad[], void *data);
// typedef void callback(double x[], int n, double f, size_t iter, void *data);

// typedef enum cg_result {tol_achieved = 0, stop_maxnfeval = 1, stop_maxiter = 2} cg_result;
// #include "cblas.h"

// int minimize_nonneg_cg(double x[restrict], int n, double *fun_val,
// 	fun_eval *obj_fun, grad_eval *grad_fun, callback *cb, void *data,
// 	double tol, size_t maxnfeval, size_t maxiter, size_t *niter, size_t *nfeval,
// 	double decr_lnsrch, double lnsrch_const, size_t max_ls,
// 	int extra_nonneg_tol, double *buffer_arr, int nthreads, int verbose)
// {
// 	double max_step;
// 	double direction_norm_sq;
// 	double grad_prev_norm_sq;
// 	double prod_grad_dir;
// 	double theta;
// 	double beta;
// 	double curr_fun_val;
// 	double new_fun_val;
// 	obj_fun(x, n, &curr_fun_val, data);
// 	*nfeval = 1;
// 	int dealloc_buffer = 0;
// 	int revert_x = 0;
// 	size_t ls;
// 	cg_result return_value = stop_maxiter;
// 	if ( maxiter <= 0 ) { maxiter = INT_MAX;}
// 	if ( maxnfeval <= 0 ) { maxnfeval = INT_MAX;}

// 	#if defined(_OPENMP) && (_OPENMP < 20080101) /* OpenMP < 3.0 */
// 	long i;
// 	long n_szt = n;
// 	#else
// 	size_t n_szt = (size_t) n;
// 	#endif

// 	/*	algorithm requires current and previous gradient and search direction, so the index
// 		at which they are written in the array is rotated each iteration to avoid unnecessary copies */
// 	int ix_rotation = 0;
// 	if (buffer_arr == NULL)
// 	{
// 		buffer_arr = (double*) malloc(sizeof(double) * n * 4);
// 		dealloc_buffer = 1;
// 		if (buffer_arr == NULL)
// 		{
// 			fprintf(stderr, "Could not allocate memory for optimization procedure\n");
// 			EXIT_FAILURE;
// 		}
// 	}
// 	double *grad_curr_n_prev = buffer_arr;
// 	double *direction_curr_n_prev = buffer_arr + 2 * n;
// 	double *restrict direction_curr = direction_curr_n_prev;
// 	double *restrict grad_curr = grad_curr_n_prev;
// 	double *restrict direction_prev;
// 	double *restrict grad_prev;

// 	/* set number of BLAS threads */
// 	#if defined(mkl_set_num_threads_local) || defined(HAS_MKL)
// 		int ignore = mkl_set_num_threads_local(nthreads);
// 	#elif defined(openblas_set_num_threads) || defined(HAS_OPENBLAS)
// 		openblas_set_num_threads(nthreads);
// 	#elif defined(_OPENMP)
// 		omp_set_num_threads(nthreads);
// 	#endif


// 	if (verbose)
// 	{
// 		printf("********************************************\n");
// 		printf("Non-negative Conjugate Gradient Optimization\n\n");
// 		printf("Number of variables to optimize: %d\n", n);
// 		if (maxiter == INT_MAX && maxnfeval == INT_MAX) {printf("[Warning: no limit on iterations and function evaluations passed]");}
// 		printf("Initial function value: %10.4f\n\n", curr_fun_val);
// 	}

// 	for (*niter = 0; *niter < maxiter; (*niter)++)
// 	{
// 		/* get gradient */
// 		grad_fun(x, n, grad_curr, data);

// 		/* determine search direction - this requires 3 passess over 'x' */

// 		/* first pass: get a capped gradient */
// 		#pragma omp parallel for schedule(static, n/nthreads) firstprivate(x, direction_curr, grad_curr, n_szt) num_threads(nthreads)
// 		for (size_t_for i = 0; i < n_szt; i++)
// 		{
// 			direction_curr[i] = (x[i] <= 0 && grad_curr[i] >= 0)? 0 : -grad_curr[i];
// 		}

// 		/* at first iteration, stop with that */
// 		if (*niter > 0)
// 		{
// 			/* second pass: calculate beta and theta constants */
// 			theta = 0;
// 			beta = 0;
// 			#pragma omp parallel for schedule(static, n/nthreads) firstprivate(x, direction_prev, grad_curr, grad_prev, n_szt) reduction(+:theta, beta) num_threads(nthreads)
// 			for (size_t_for i = 0; i < n_szt; i++)
// 			{
// 				theta += ( x[i] <= 0 )? 0 : grad_curr[i] * direction_prev[i];
// 				beta += ( x[i] <= 0 )? 0 : grad_curr[i] * (grad_curr[i] - grad_prev[i]);
// 			}
// 			theta /= grad_prev_norm_sq;
// 			beta /= grad_prev_norm_sq;

// 			/* third pass: add to direction info on previous direction and gradient differences */
// 			#pragma omp parallel for schedule(static, n/nthreads) firstprivate(x, direction_curr, direction_prev, grad_curr, grad_prev, n_szt, theta, beta) num_threads(nthreads)
// 			for (size_t_for i = 0; i < n_szt; i++)
// 			{
// 				direction_curr[i] += ( x[i] <= 0 )? 0 : beta * direction_prev[i] - theta * (grad_curr[i] - grad_prev[i]);
// 			}

// 		}

// 		/* check if stop criterion is satisfied */
// 		prod_grad_dir = cblas_ddot(n, grad_curr, 1, direction_curr, 1);
// 		if ( fabs(prod_grad_dir) <= tol )
// 		{
// 			return_value = tol_achieved;
// 			goto terminate_procedure;
// 		}

// 		/* determine maximum step size */
// 		max_step = 1;
// 		#if defined(_OPENMP)
// 		#pragma omp parallel for schedule(static, n/nthreads) firstprivate(x, direction_curr, n_szt) reduction(min: max_step) num_threads(nthreads)
// 		for (size_t_for i = 0; i < n_szt; i++)
// 		{
// 			max_step = (direction_curr[i] < 0)? -x[i] / direction_curr[i] : 1;
// 		}
// 		max_step = min(max_step, 1);

// 		#else
// 		for (size_t i = 0; i < n_szt; i++)
// 		{
// 			if (direction_curr[i] < 0) { max_step = min(max_step, -x[i] / direction_curr[i]); }
// 		}
// 		#endif

// 		/* perform line search */
// 		cblas_daxpy(n, max_step, direction_curr, 1, x, 1);
// 		direction_norm_sq = cblas_ddot(n, direction_curr, 1, direction_curr, 1);
// 		if (extra_nonneg_tol)
// 		{
// 			#pragma omp parallel for schedule(static, n/nthreads) firstprivate(x, n_szt) num_threads(nthreads)
// 			for (size_t_for i = 0; i < n_szt; i++){x[i] = (x[i] <= 0)? 0 : x[i];}
// 		}
// 		for (ls = 0; ls < max_ls; ls++)
// 		{
// 			obj_fun(x, n, &new_fun_val, data);
// 			if ( !isinf(new_fun_val) && !isnan(new_fun_val) )
// 			{
// 				if (new_fun_val <=  curr_fun_val - lnsrch_const * square(max_step * pow(decr_lnsrch, ls)) * direction_norm_sq)
// 					{ break; }
// 			}
// 			(*nfeval)++; if (*nfeval >= maxnfeval) { revert_x = 1; return_value = stop_maxnfeval; goto terminate_procedure; }
// 			/* go to new step size by modifying x in-place */
// 			cblas_daxpy(n, max_step * ( pow(decr_lnsrch, ls + 1) - pow(decr_lnsrch, ls) ), direction_curr, 1, x, 1);
// 		}
// 		curr_fun_val = new_fun_val;
// 		if ( cb != NULL) { cb(x, n, curr_fun_val, *niter, data); }

// 		/* update norm of gradient */
// 		grad_prev_norm_sq = cblas_ddot(n, grad_curr, 1, grad_curr, 1);

// 		/* next time, write to the other side of grad and dir arrays */
// 		direction_prev = direction_curr;
// 		grad_prev = grad_curr;
// 		ix_rotation = incr_ix_rotation(ix_rotation);
// 		direction_curr = direction_curr_n_prev + get_curr_ix_rotation(ix_rotation, n);
// 		grad_curr = grad_curr_n_prev + get_curr_ix_rotation(ix_rotation, n);
// 		if (verbose)
// 		{
// 			printf("Iteration %3d : f(x) = %10.4f, |<g(x), d(x)>| = %12.4f, nfev = %3d, ls = %2d\n",
// 					(int) *niter + 1, curr_fun_val, fabs(prod_grad_dir), (int) *nfeval, (int) ls +1 );
// 		}
// 	}

// 	terminate_procedure:
// 		if (dealloc_buffer) { free(buffer_arr); }
// 		if (revert_x) { cblas_daxpy(n, -max_step * pow(decr_lnsrch, ls), direction_curr, 1, x, 1); }
// 		if (verbose)
// 		{
// 			if (return_value == tol_achieved) 	{ printf("\nTerminated: |<g(x), d(x)>| driven below tol.\n"); }
// 			if (return_value == stop_maxnfeval) { printf("\nTerminated: reached maximum number of function evaluations\n"); }
// 			if (return_value == stop_maxiter) 	{ printf("\nTerminated: reached maximum number of iterations\n"); }
// 			printf("Last f(x) = %10.4f\n\n", curr_fun_val);
// 		}
// 		*fun_val = curr_fun_val;
// 	return (int) return_value;
// }





/*	Example usage of non-negative conjugate gradient optimizer
	with Rosenbrock's (a.k.a. banana) function.
	Implementation taken from here:
	https://gist.github.com/Bismarrck/772d4ae4a8c5ddcd37b8
*/
#include <stdlib.h>
#include <stdio.h>
// #include "nonnegcg.h"
typedef void fun_eval(double x[], int n, double *f, void *data);
typedef void grad_eval(double x[], int n, double grad[], void *data);
typedef void callback(double x[], int n, double f, size_t iter, void *data);

typedef enum cg_result {tol_achieved = 0, stop_maxnfeval = 1, stop_maxiter = 2} cg_result;

int minimize_nonneg_cg(double x[restrict], int n, double *fun_val,
	fun_eval *obj_fun, grad_eval *grad_fun, callback *cb, void *data,
	double tol, size_t maxnfeval, size_t maxiter, size_t *niter, size_t *nfeval,
	double decr_lnsrch, double lnsrch_const, size_t max_ls,
	int extra_nonneg_tol, double *buffer_arr, int nthreads, int verbose);

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

	// minimize_nonneg_cg(x0, n, &fun_val,
	// 	rosen, rosen_der, NULL, NULL,
	// 	tol, maxnfeval, maxiter, &niter, &nfeval,
	// 	decr_lnsrch, lnsrch_const, max_ls,
	// 	extra_nonneg_tol, NULL, nthreads, verbose);
	double *buffer = (double*) malloc(sizeof(double) * 4 * 4);

	minimize_nonneg_cg(x, n, &fun_val,
		rosen, rosen_der, NULL, NULL,
		tol, maxnfeval, maxiter, &niter, &nfeval,
		decr_lnsrch, lnsrch_const, max_ls,
		extra_nonneg_tol, buffer, nthreads, verbose);

	printf("Optimal values of x: [ ");
	for (size_t i = 0; i < n; i++){printf("%f ", x[i]);} printf("]\n");

	return 0;
}
