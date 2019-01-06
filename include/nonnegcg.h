#ifdef __cplusplus
extern "C" {
#endif

/*	Functions passed to the optimizer must have these propotypes */
typedef void fun_eval(double x[], int n, double *f, void *data);
typedef void grad_eval(double x[], int n, double grad[], void *data);
typedef void callback(double x[], int n, double f, size_t iter, void *data);

/*	Optimizer return codes - note that they are return as int and not as enum */
typedef enum cg_result {tol_achieved = 0, stop_maxnfeval = 1, stop_maxiter = 2} cg_result;

/*	Non-negative conjugate gradient optimizer
	
	Minimizes a function subject to non-negativity constraints on all the variables,
	using a modified Polak-Rubiere-Polyak conjugate gradient method. Implementation
	is based on the paper:

	Li, C. (2013). A conjugate gradient type method for the nonnegative constraints optimization problems. Journal of Applied Mathematics, 2013.
	
	x (in, out)		: At input, starting point (must be a feasible point). At output, optimal values calculated by the optimizer.
	n 				: Number of variables in the optimization problem
	fun_val (out)	: Value of the function achieved at the end of the procedure
	obj_fun			: function that calculates the objective value (must be written into the *f pointer passed to it)
	grad_fun		: function that calculates the gradient (must be written into the grad[] array passed to it)
	cb				: callback function to execute at the end of each iteration
	data			: Extra data to pass to the functions that evaluate objective, gradient, and callback (must be cast to void pointer)
	tol				: Tolerance for <gradient, direction>
					  (Recommended: <1e-3)
	maxnfeval		: Maximum number of function evaluations
					  (Recommended: >1000)
	maxiter			: Maximum number of CG iterations to run
					  (Recommended: >100, but note that steps are always feasible descent directions)
	niter (out)		: Number of CG iterations performed
	nfeval (out)	: Number of function evaluations performed
	decr_lnsrch		: Number by which to decrease the step size after each unsuccessful line search
					  (Recommended: 0.5)
	lnsrch_const	: Acceptance parameter for the line search procedure
					  (Recommended: 0.01)
	max_ls			: Maximum number of line search trials per iteration
					  (Recommended: 20)
	extra_nonneg_tol: Ensure extra non-negative tolerance by explicitly setting elements that are <=0 to zero at each iteration
					  (Recommended: 0)
	buffer_arr		: Array of dimensions (4*n). Will allocate it and then free it if passing NULL.
	nthreads		: Number of parallel threads to use
	verbose			: Whether to print convergence messages
*/
int minimize_nonneg_cg(double x[], int n, double *fun_val,
	fun_eval *obj_fun, grad_eval *grad_fun, callback *cb, void *data,
	double tol, size_t maxnfeval, size_t maxiter, size_t *niter, size_t *nfeval,
	double decr_lnsrch, double lnsrch_const, size_t max_ls,
	int extra_nonneg_tol, double *buffer_arr, int nthreads, int verbose);

#ifdef __cplusplus
}
#endif
