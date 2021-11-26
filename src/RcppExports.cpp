// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cal2norm
double cal2norm(const arma::mat A, const arma::umat gIdx, const bool Atranspose, const int TgCB_T);
RcppExport SEXP _mixlasso_cal2norm(SEXP ASEXP, SEXP gIdxSEXP, SEXP AtransposeSEXP, SEXP TgCB_TSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::umat >::type gIdx(gIdxSEXP);
    Rcpp::traits::input_parameter< const bool >::type Atranspose(AtransposeSEXP);
    Rcpp::traits::input_parameter< const int >::type TgCB_T(TgCB_TSEXP);
    rcpp_result_gen = Rcpp::wrap(cal2norm(A, gIdx, Atranspose, TgCB_T));
    return rcpp_result_gen;
END_RCPP
}
// shrink
arma::mat shrink(const arma::mat A, const arma::umat gIdx, const bool Atranspose, const int TgCB_T);
RcppExport SEXP _mixlasso_shrink(SEXP ASEXP, SEXP gIdxSEXP, SEXP AtransposeSEXP, SEXP TgCB_TSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::umat >::type gIdx(gIdxSEXP);
    Rcpp::traits::input_parameter< const bool >::type Atranspose(AtransposeSEXP);
    Rcpp::traits::input_parameter< const int >::type TgCB_T(TgCB_TSEXP);
    rcpp_result_gen = Rcpp::wrap(shrink(A, gIdx, Atranspose, TgCB_T));
    return rcpp_result_gen;
END_RCPP
}
// XXeigen
arma::vec XXeigen(const arma::mat& X, const arma::mat& V, const arma::uvec& tIdx);
RcppExport SEXP _mixlasso_XXeigen(SEXP XSEXP, SEXP VSEXP, SEXP tIdxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type tIdx(tIdxSEXP);
    rcpp_result_gen = Rcpp::wrap(XXeigen(X, V, tIdx));
    return rcpp_result_gen;
END_RCPP
}
// mixLassoLoop
arma::sp_mat mixLassoLoop(const arma::mat& V, const arma::sp_mat& C, const arma::umat& gIdx, const bool Tglasso, const arma::uvec& tIdx, const arma::mat& X, const arma::mat& Y, const bool intercept, const int num_nonpen, arma::vec L, const double lambda, const arma::vec& option, const double mu, const int NoVar, const double gamma, const arma::umat& y_mis, const double alpha);
RcppExport SEXP _mixlasso_mixLassoLoop(SEXP VSEXP, SEXP CSEXP, SEXP gIdxSEXP, SEXP TglassoSEXP, SEXP tIdxSEXP, SEXP XSEXP, SEXP YSEXP, SEXP interceptSEXP, SEXP num_nonpenSEXP, SEXP LSEXP, SEXP lambdaSEXP, SEXP optionSEXP, SEXP muSEXP, SEXP NoVarSEXP, SEXP gammaSEXP, SEXP y_misSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type gIdx(gIdxSEXP);
    Rcpp::traits::input_parameter< const bool >::type Tglasso(TglassoSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type tIdx(tIdxSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const bool >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< const int >::type num_nonpen(num_nonpenSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type L(LSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type option(optionSEXP);
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const int >::type NoVar(NoVarSEXP);
    Rcpp::traits::input_parameter< const double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type y_mis(y_misSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(mixLassoLoop(V, C, gIdx, Tglasso, tIdx, X, Y, intercept, num_nonpen, L, lambda, option, mu, NoVar, gamma, y_mis, alpha));
    return rcpp_result_gen;
END_RCPP
}
// treeLassoLoop
arma::sp_mat treeLassoLoop(const arma::mat& X, const arma::mat& Y, const arma::sp_mat& C, const arma::umat& gIdx, const double TauNorm, const bool intercept, const int num_nonpen, const double lambda, const arma::vec& option, const double mu);
RcppExport SEXP _mixlasso_treeLassoLoop(SEXP XSEXP, SEXP YSEXP, SEXP CSEXP, SEXP gIdxSEXP, SEXP TauNormSEXP, SEXP interceptSEXP, SEXP num_nonpenSEXP, SEXP lambdaSEXP, SEXP optionSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type gIdx(gIdxSEXP);
    Rcpp::traits::input_parameter< const double >::type TauNorm(TauNormSEXP);
    Rcpp::traits::input_parameter< const bool >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< const int >::type num_nonpen(num_nonpenSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type option(optionSEXP);
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(treeLassoLoop(X, Y, C, gIdx, TauNorm, intercept, num_nonpen, lambda, option, mu));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mixlasso_cal2norm", (DL_FUNC) &_mixlasso_cal2norm, 4},
    {"_mixlasso_shrink", (DL_FUNC) &_mixlasso_shrink, 4},
    {"_mixlasso_XXeigen", (DL_FUNC) &_mixlasso_XXeigen, 3},
    {"_mixlasso_mixLassoLoop", (DL_FUNC) &_mixlasso_mixLassoLoop, 17},
    {"_mixlasso_treeLassoLoop", (DL_FUNC) &_mixlasso_treeLassoLoop, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_mixlasso(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
