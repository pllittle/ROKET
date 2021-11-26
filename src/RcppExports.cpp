// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/ROKET.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Rcpp_run_OT
arma::mat Rcpp_run_OT(const arma::vec& XX, const arma::vec& YY, const arma::mat& COST_XY, const double& EPS, const double& LAMBDA1, const double& LAMBDA2, const bool& balance, const bool& highLAM_lowMU, const double& conv, const arma::uword& max_iter, const bool& show, const arma::uword& show_iter);
static SEXP _ROKET_Rcpp_run_OT_try(SEXP XXSEXP, SEXP YYSEXP, SEXP COST_XYSEXP, SEXP EPSSEXP, SEXP LAMBDA1SEXP, SEXP LAMBDA2SEXP, SEXP balanceSEXP, SEXP highLAM_lowMUSEXP, SEXP convSEXP, SEXP max_iterSEXP, SEXP showSEXP, SEXP show_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
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
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _ROKET_Rcpp_run_OT(SEXP XXSEXP, SEXP YYSEXP, SEXP COST_XYSEXP, SEXP EPSSEXP, SEXP LAMBDA1SEXP, SEXP LAMBDA2SEXP, SEXP balanceSEXP, SEXP highLAM_lowMUSEXP, SEXP convSEXP, SEXP max_iterSEXP, SEXP showSEXP, SEXP show_iterSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_ROKET_Rcpp_run_OT_try(XXSEXP, YYSEXP, COST_XYSEXP, EPSSEXP, LAMBDA1SEXP, LAMBDA2SEXP, balanceSEXP, highLAM_lowMUSEXP, convSEXP, max_iterSEXP, showSEXP, show_iterSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// Rcpp_run_full_OT
Rcpp::List Rcpp_run_full_OT(const arma::mat& COST, const arma::mat& ZZ, const double& EPS, const double& LAMBDA1, const double& LAMBDA2, const bool& balance, const bool& highLAM_lowMU, const double& conv, const arma::uword& max_iter, const int& ncores, const bool& show, const arma::uword& show_iter);
static SEXP _ROKET_Rcpp_run_full_OT_try(SEXP COSTSEXP, SEXP ZZSEXP, SEXP EPSSEXP, SEXP LAMBDA1SEXP, SEXP LAMBDA2SEXP, SEXP balanceSEXP, SEXP highLAM_lowMUSEXP, SEXP convSEXP, SEXP max_iterSEXP, SEXP ncoresSEXP, SEXP showSEXP, SEXP show_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
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
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _ROKET_Rcpp_run_full_OT(SEXP COSTSEXP, SEXP ZZSEXP, SEXP EPSSEXP, SEXP LAMBDA1SEXP, SEXP LAMBDA2SEXP, SEXP balanceSEXP, SEXP highLAM_lowMUSEXP, SEXP convSEXP, SEXP max_iterSEXP, SEXP ncoresSEXP, SEXP showSEXP, SEXP show_iterSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_ROKET_Rcpp_run_full_OT_try(COSTSEXP, ZZSEXP, EPSSEXP, LAMBDA1SEXP, LAMBDA2SEXP, balanceSEXP, highLAM_lowMUSEXP, convSEXP, max_iterSEXP, ncoresSEXP, showSEXP, show_iterSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// Rcpp_KernTest
Rcpp::List Rcpp_KernTest(const arma::vec& RESI, const Rcpp::List& KK, const arma::uword& nPERMS, const arma::uword& iter1, const arma::uword& iter2, const bool& verbose);
static SEXP _ROKET_Rcpp_KernTest_try(SEXP RESISEXP, SEXP KKSEXP, SEXP nPERMSSEXP, SEXP iter1SEXP, SEXP iter2SEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type RESI(RESISEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type KK(KKSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type nPERMS(nPERMSSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type iter1(iter1SEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type iter2(iter2SEXP);
    Rcpp::traits::input_parameter< const bool& >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_KernTest(RESI, KK, nPERMS, iter1, iter2, verbose));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _ROKET_Rcpp_KernTest(SEXP RESISEXP, SEXP KKSEXP, SEXP nPERMSSEXP, SEXP iter1SEXP, SEXP iter2SEXP, SEXP verboseSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_ROKET_Rcpp_KernTest_try(RESISEXP, KKSEXP, nPERMSSEXP, iter1SEXP, iter2SEXP, verboseSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _ROKET_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("arma::mat(*Rcpp_run_OT)(const arma::vec&,const arma::vec&,const arma::mat&,const double&,const double&,const double&,const bool&,const bool&,const double&,const arma::uword&,const bool&,const arma::uword&)");
        signatures.insert("Rcpp::List(*Rcpp_run_full_OT)(const arma::mat&,const arma::mat&,const double&,const double&,const double&,const bool&,const bool&,const double&,const arma::uword&,const int&,const bool&,const arma::uword&)");
        signatures.insert("Rcpp::List(*Rcpp_KernTest)(const arma::vec&,const Rcpp::List&,const arma::uword&,const arma::uword&,const arma::uword&,const bool&)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _ROKET_RcppExport_registerCCallable() { 
    R_RegisterCCallable("ROKET", "_ROKET_Rcpp_run_OT", (DL_FUNC)_ROKET_Rcpp_run_OT_try);
    R_RegisterCCallable("ROKET", "_ROKET_Rcpp_run_full_OT", (DL_FUNC)_ROKET_Rcpp_run_full_OT_try);
    R_RegisterCCallable("ROKET", "_ROKET_Rcpp_KernTest", (DL_FUNC)_ROKET_Rcpp_KernTest_try);
    R_RegisterCCallable("ROKET", "_ROKET_RcppExport_validate", (DL_FUNC)_ROKET_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_ROKET_Rcpp_run_OT", (DL_FUNC) &_ROKET_Rcpp_run_OT, 12},
    {"_ROKET_Rcpp_run_full_OT", (DL_FUNC) &_ROKET_Rcpp_run_full_OT, 12},
    {"_ROKET_Rcpp_KernTest", (DL_FUNC) &_ROKET_Rcpp_KernTest, 6},
    {"_ROKET_RcppExport_registerCCallable", (DL_FUNC) &_ROKET_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_ROKET(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
