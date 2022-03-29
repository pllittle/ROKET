// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Rcpp_run_OT
arma::mat Rcpp_run_OT(const arma::vec& XX, const arma::vec& YY, const arma::mat& COST_XY, const double& EPS, const double& LAMBDA1, const double& LAMBDA2, const bool& balance, const bool& highLAM_lowMU, const double& conv, const arma::uword& max_iter, const bool& show, const arma::uword& show_iter);
RcppExport SEXP _ROKET_Rcpp_run_OT(SEXP XXSEXP, SEXP YYSEXP, SEXP COST_XYSEXP, SEXP EPSSEXP, SEXP LAMBDA1SEXP, SEXP LAMBDA2SEXP, SEXP balanceSEXP, SEXP highLAM_lowMUSEXP, SEXP convSEXP, SEXP max_iterSEXP, SEXP showSEXP, SEXP show_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type XX(XXSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type YY(YYSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type COST_XY(COST_XYSEXP);
    Rcpp::traits::input_parameter< const double& >::type EPS(EPSSEXP);
    Rcpp::traits::input_parameter< const double& >::type LAMBDA1(LAMBDA1SEXP);
    Rcpp::traits::input_parameter< const double& >::type LAMBDA2(LAMBDA2SEXP);
    Rcpp::traits::input_parameter< const bool& >::type balance(balanceSEXP);
    Rcpp::traits::input_parameter< const bool& >::type highLAM_lowMU(highLAM_lowMUSEXP);
    Rcpp::traits::input_parameter< const double& >::type conv(convSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const bool& >::type show(showSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type show_iter(show_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_run_OT(XX, YY, COST_XY, EPS, LAMBDA1, LAMBDA2, balance, highLAM_lowMU, conv, max_iter, show, show_iter));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_run_full_OT
Rcpp::List Rcpp_run_full_OT(const arma::mat& COST, const arma::mat& ZZ, const double& EPS, const double& LAMBDA1, const double& LAMBDA2, const bool& balance, const bool& highLAM_lowMU, const double& conv, const arma::uword& max_iter, const int& ncores, const bool& show, const arma::uword& show_iter);
RcppExport SEXP _ROKET_Rcpp_run_full_OT(SEXP COSTSEXP, SEXP ZZSEXP, SEXP EPSSEXP, SEXP LAMBDA1SEXP, SEXP LAMBDA2SEXP, SEXP balanceSEXP, SEXP highLAM_lowMUSEXP, SEXP convSEXP, SEXP max_iterSEXP, SEXP ncoresSEXP, SEXP showSEXP, SEXP show_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type COST(COSTSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ZZ(ZZSEXP);
    Rcpp::traits::input_parameter< const double& >::type EPS(EPSSEXP);
    Rcpp::traits::input_parameter< const double& >::type LAMBDA1(LAMBDA1SEXP);
    Rcpp::traits::input_parameter< const double& >::type LAMBDA2(LAMBDA2SEXP);
    Rcpp::traits::input_parameter< const bool& >::type balance(balanceSEXP);
    Rcpp::traits::input_parameter< const bool& >::type highLAM_lowMU(highLAM_lowMUSEXP);
    Rcpp::traits::input_parameter< const double& >::type conv(convSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< const bool& >::type show(showSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type show_iter(show_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_run_full_OT(COST, ZZ, EPS, LAMBDA1, LAMBDA2, balance, highLAM_lowMU, conv, max_iter, ncores, show, show_iter));
    return rcpp_result_gen;
END_RCPP
}
// Rcpp_KernTest
Rcpp::List Rcpp_KernTest(const arma::vec& RESI, const arma::cube& cKK, const arma::umat& OMNI, const arma::uword& nPERMS);
RcppExport SEXP _ROKET_Rcpp_KernTest(SEXP RESISEXP, SEXP cKKSEXP, SEXP OMNISEXP, SEXP nPERMSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type RESI(RESISEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type cKK(cKKSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type OMNI(OMNISEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type nPERMS(nPERMSSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_KernTest(RESI, cKK, OMNI, nPERMS));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ROKET_Rcpp_run_OT", (DL_FUNC) &_ROKET_Rcpp_run_OT, 12},
    {"_ROKET_Rcpp_run_full_OT", (DL_FUNC) &_ROKET_Rcpp_run_full_OT, 12},
    {"_ROKET_Rcpp_KernTest", (DL_FUNC) &_ROKET_Rcpp_KernTest, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_ROKET(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
