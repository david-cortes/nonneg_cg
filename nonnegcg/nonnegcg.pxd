ctypedef void fun_eval(double *x, int n, double *f, void *data)
ctypedef void grad_eval(double *x, int n, double *grad, void *data)
ctypedef void callback(double *x, int n, size_t iter, void *data) 
ctypedef enum cg_result:
		tol_achieved = 0
		stop_maxnfeval = 1
		stop_maxiter = 2
int minimize_nonneg_cg(double *x, int n, double *fun_val,
	fun_eval *obj_fun, grad_eval *grad_fun, callback *cb, void *data,
	double tol, size_t maxnfeval, size_t maxiter, size_t *niter, size_t *nfeval,
	double decr_lnsrch, double lnsrch_const, size_t max_ls,
	int extra_nonneg_tol, double *buffer_arr, int nthreads, int verbose)
