// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// run_minimizer
Rcpp::List run_minimizer(Rcpp::Function feval, Rcpp::Function geval, Rcpp::Function docall, Rcpp::NumericVector x0, Rcpp::List fargs, double tol, size_t maxnfeval, size_t maxiter, double decr_lnsrch, double lnsrch_const, size_t max_ls, int extra_nonneg_tol, int nthreads, int verbose);
RcppExport SEXP _nonneg_cg_run_minimizer(SEXP fevalSEXP, SEXP gevalSEXP, SEXP docallSEXP, SEXP x0SEXP, SEXP fargsSEXP, SEXP tolSEXP, SEXP maxnfevalSEXP, SEXP maxiterSEXP, SEXP decr_lnsrchSEXP, SEXP lnsrch_constSEXP, SEXP max_lsSEXP, SEXP extra_nonneg_tolSEXP, SEXP nthreadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Function >::type feval(fevalSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type geval(gevalSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type docall(docallSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type fargs(fargsSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< size_t >::type maxnfeval(maxnfevalSEXP);
    Rcpp::traits::input_parameter< size_t >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type decr_lnsrch(decr_lnsrchSEXP);
    Rcpp::traits::input_parameter< double >::type lnsrch_const(lnsrch_constSEXP);
    Rcpp::traits::input_parameter< size_t >::type max_ls(max_lsSEXP);
    Rcpp::traits::input_parameter< int >::type extra_nonneg_tol(extra_nonneg_tolSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(run_minimizer(feval, geval, docall, x0, fargs, tol, maxnfeval, maxiter, decr_lnsrch, lnsrch_const, max_ls, extra_nonneg_tol, nthreads, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_nonneg_cg_run_minimizer", (DL_FUNC) &_nonneg_cg_run_minimizer, 14},
    {NULL, NULL, 0}
};

void R_init_nonnegcg2(DllInfo *info);
RcppExport void R_init_nonneg_cg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_init_nonnegcg2(dll);
}
