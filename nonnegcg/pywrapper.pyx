import numpy as np
cimport numpy as np
from libc.string cimport memcpy
from scipy.optimize import OptimizeResult
import multiprocessing

cdef extern from "../src/nonnegcg.c":
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

cdef void wrapped_eval_f(double *x, int n, double *f, void *data):
	py_tuple = <tuple> data
	f[0] = py_tuple[0](py_tuple[2], *py_tuple[5])

cdef void wrapped_eval_g(double *x, int n, double *grad, void *data):
	py_tuple = <tuple> data
	if py_tuple[4][0] == 0:
		py_tuple[3][: n] = py_tuple[1](py_tuple[2], *py_tuple[5])
		py_tuple[4][0] = 1
	else:
		py_tuple[3][n : 2*n] = py_tuple[1](py_tuple[2], *py_tuple[5])
		py_tuple[4][0] = 0

def minimize_nncg(x0, fun, grad, args=(), decr_lnsrch=.5, lnsrch_const=.01, maxiter=200,
	maxnfeval=1500, max_ls=20, tol=1e-4, extra_nonneg_tol=False, nthreads=-1, verbose=False):
	"""
	Minimize a function subject to non-negativity constraints on all variables

	Note
	----
	The C function can also be used directly in Cython by cimport'ing it (see .pxd header), and in C
	by linking  a Python extension to `nonnegcg.minimize_nncg.__file__` and including `.h` header that
	this package will also intall.

	Parameters
	----------
	x0 : array(n_variables,)
		Starting point. Must be a feasible point (>=0). Be aware that it might be modified in-place.
	fun : function(x, *args) -> float
		Function that evaluates the function to minimize at 'x'
	grad : function(x, *args) -> array(n_variables,)
		Function that evaluates the gradient of the function to minimze at 'x'
	args : tuple
		Extra arguments to pass to 'fun' and 'grad'
	decr_lnsrch : float(0,1)
		Number by which to decrease the step size after each unsuccessful line search
	lnsrch_const : float(0,1)
		Acceptance parameter for the line search procedure
	maxiter : int > 0
		Maximum number of CG iterations to run
	maxnfeval : int > 0
		Maximum number of function evaluations
	max_ls : int > 0
		Maximum number of line search trials per iteration
	tol : float > 0
		Tolerance for <gradient, direction>
	extra_nonneg_tol : bool
		Ensure extra non-negative tolerance by explicitly setting elements that are <=0 to zero at each iteration
	nthreads : int
		Number of parallel threads to use
	verbose : bool
		Whether to print convergence messages.
		Note that if you are running this in ipython, the messages will appear in the console but not on the notebook.

	Returns
	-------
	result : OptimizeResult
		A dict-like structure with the final variables and info on the optimization procedure

	References
	----------
	Li, C. (2013). A conjugate gradient type method for the nonnegative constraints optimization problems. Journal of Applied Mathematics, 2013.
	"""
	assert isinstance(args, tuple)
	if nthreads < 1:
		nthreads = multiprocessing.cpu_count()
		if nthreads is None:
			nthreads = 1 
	assert nthreads > 0
	assert isinstance(nthreads, int)
	assert tol > 0
	assert (decr_lnsrch > 0) and (decr_lnsrch < 1)
	assert (lnsrch_const > 0) and (lnsrch_const < 1)
	assert maxiter > 0
	assert maxnfeval > 0
	assert max_ls > 0

	cdef np.ndarray[double, ndim=1] x_np = np.array(x0).astype('float64').reshape(-1)
	if x_np.min() < 0:
		raise ValueError("'x0' must be a feasible point.")
	cdef double c_decr_lnsrch = <double> decr_lnsrch
	cdef double c_lnsrch_const = <double> lnsrch_const
	cdef double c_tol = <double> tol
	cdef size_t c_maxiter = <size_t> maxiter
	cdef size_t c_maxnfeval = <size_t> maxnfeval
	cdef size_t c_max_ls = <size_t> max_ls
	cdef int c_extra_nonneg_tol = <int> bool(extra_nonneg_tol)
	cdef int c_verbose = <int> bool(verbose)
	cdef int c_nthreads = <int> nthreads
	cdef int n = x_np.shape[0]

	cdef double slot_f
	cdef size_t slot_niter, slot_nfeval
	cdef np.ndarray[double, ndim=1] slot_buffer = np.zeros(n * 4, dtype='float64')
	cdef np.ndarray[double, ndim=1] grad_cycler = np.zeros(1, dtype='float64')
	tuple_give = (fun, grad, x_np, slot_buffer, grad_cycler, args)

	cdef cg_result res = <cg_result> minimize_nonneg_cg(&x_np[0], n, &slot_f,
		wrapped_eval_f, wrapped_eval_g, NULL, <void*> tuple_give,
		c_tol, c_maxnfeval, c_maxiter, &slot_niter, &slot_nfeval,
		c_decr_lnsrch, c_lnsrch_const, c_max_ls,
		c_extra_nonneg_tol, &slot_buffer[0], c_nthreads, c_verbose)

	dct_out = dict()
	dct_out['fun'] = slot_f
	dct_out['status'] = {0:'tol_achieved', 1:'stop_maxnfeval', 2:'stop_maxiter'}[res]
	dct_out['nfev'] = slot_nfeval
	dct_out['nit'] = slot_niter
	dct_out['x'] = x_np

	return OptimizeResult(dct_out)
