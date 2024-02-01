// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gpxi_int
double gpxi_int(int i, int m, NumericVector beta, NumericVector S0, NumericVector S1, double sig2);
RcppExport SEXP _SurrogateBMA_gpxi_int(SEXP iSEXP, SEXP mSEXP, SEXP betaSEXP, SEXP S0SEXP, SEXP S1SEXP, SEXP sig2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type S1(S1SEXP);
    Rcpp::traits::input_parameter< double >::type sig2(sig2SEXP);
    rcpp_result_gen = Rcpp::wrap(gpxi_int(i, m, beta, S0, S1, sig2));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _SurrogateBMA_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SurrogateBMA_gpxi_int", (DL_FUNC) &_SurrogateBMA_gpxi_int, 6},
    {"_SurrogateBMA_rcpp_hello_world", (DL_FUNC) &_SurrogateBMA_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_SurrogateBMA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}