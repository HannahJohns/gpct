// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getRsdt
NumericMatrix getRsdt(const IntegerMatrix& pairwise, const NumericVector& group, const NumericVector& weight, const double tol);
RcppExport SEXP _gpct_getRsdt(SEXP pairwiseSEXP, SEXP groupSEXP, SEXP weightSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type pairwise(pairwiseSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type group(groupSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(getRsdt(pairwise, group, weight, tol));
    return rcpp_result_gen;
END_RCPP
}
// get_nullSE
List get_nullSE(const NumericMatrix& winProb_X, const NumericMatrix& winProb_Y, const NumericVector& weight_X, const NumericVector& weight_Y);
RcppExport SEXP _gpct_get_nullSE(SEXP winProb_XSEXP, SEXP winProb_YSEXP, SEXP weight_XSEXP, SEXP weight_YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type winProb_X(winProb_XSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type winProb_Y(winProb_YSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weight_X(weight_XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weight_Y(weight_YSEXP);
    rcpp_result_gen = Rcpp::wrap(get_nullSE(winProb_X, winProb_Y, weight_X, weight_Y));
    return rcpp_result_gen;
END_RCPP
}
// get_Rcc
NumericMatrix get_Rcc(const IntegerMatrix& pairwise, const NumericVector& group, const NumericVector& weight, const double tol);
RcppExport SEXP _gpct_get_Rcc(SEXP pairwiseSEXP, SEXP groupSEXP, SEXP weightSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type pairwise(pairwiseSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type group(groupSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(get_Rcc(pairwise, group, weight, tol));
    return rcpp_result_gen;
END_RCPP
}
// get_pw_numeric
IntegerMatrix get_pw_numeric(const NumericMatrix& outcome, const NumericVector& discriminant);
RcppExport SEXP _gpct_get_pw_numeric(SEXP outcomeSEXP, SEXP discriminantSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type outcome(outcomeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type discriminant(discriminantSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pw_numeric(outcome, discriminant));
    return rcpp_result_gen;
END_RCPP
}
// get_pw_surv
IntegerMatrix get_pw_surv(const IntegerMatrix& outcome_exists, const NumericMatrix& outcome, const NumericVector& direction, const bool censored_as_tie, const double minDiff);
RcppExport SEXP _gpct_get_pw_surv(SEXP outcome_existsSEXP, SEXP outcomeSEXP, SEXP directionSEXP, SEXP censored_as_tieSEXP, SEXP minDiffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type outcome_exists(outcome_existsSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type outcome(outcomeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type direction(directionSEXP);
    Rcpp::traits::input_parameter< const bool >::type censored_as_tie(censored_as_tieSEXP);
    Rcpp::traits::input_parameter< const double >::type minDiff(minDiffSEXP);
    rcpp_result_gen = Rcpp::wrap(get_pw_surv(outcome_exists, outcome, direction, censored_as_tie, minDiff));
    return rcpp_result_gen;
END_RCPP
}
// countInconsistency
List countInconsistency(IntegerMatrix x, NumericVector weights, bool includeTies, int m);
RcppExport SEXP _gpct_countInconsistency(SEXP xSEXP, SEXP weightsSEXP, SEXP includeTiesSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type includeTies(includeTiesSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(countInconsistency(x, weights, includeTies, m));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gpct_getRsdt", (DL_FUNC) &_gpct_getRsdt, 4},
    {"_gpct_get_nullSE", (DL_FUNC) &_gpct_get_nullSE, 4},
    {"_gpct_get_Rcc", (DL_FUNC) &_gpct_get_Rcc, 4},
    {"_gpct_get_pw_numeric", (DL_FUNC) &_gpct_get_pw_numeric, 2},
    {"_gpct_get_pw_surv", (DL_FUNC) &_gpct_get_pw_surv, 5},
    {"_gpct_countInconsistency", (DL_FUNC) &_gpct_countInconsistency, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_gpct(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
