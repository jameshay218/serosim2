// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// likelihood_titre
double likelihood_titre(NumericVector expected, NumericVector data, double error, int MAX_TITRE);
RcppExport SEXP _serosim2_likelihood_titre(SEXP expectedSEXP, SEXP dataSEXP, SEXP errorSEXP, SEXP MAX_TITRESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type expected(expectedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type error(errorSEXP);
    Rcpp::traits::input_parameter< int >::type MAX_TITRE(MAX_TITRESEXP);
    rcpp_result_gen = Rcpp::wrap(likelihood_titre(expected, data, error, MAX_TITRE));
    return rcpp_result_gen;
END_RCPP
}
// vector_sort
NumericVector vector_sort(NumericVector x);
RcppExport SEXP _serosim2_vector_sort(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(vector_sort(x));
    return rcpp_result_gen;
END_RCPP
}
// vector_order
IntegerVector vector_order(NumericVector x);
RcppExport SEXP _serosim2_vector_order(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(vector_order(x));
    return rcpp_result_gen;
END_RCPP
}
// individual_sim
NumericVector individual_sim(NumericVector mu_pars, NumericVector tp_pars, NumericVector m_pars, NumericVector ti_pars, double y0b, double lower_titre_bound, NumericVector times);
RcppExport SEXP _serosim2_individual_sim(SEXP mu_parsSEXP, SEXP tp_parsSEXP, SEXP m_parsSEXP, SEXP ti_parsSEXP, SEXP y0bSEXP, SEXP lower_titre_boundSEXP, SEXP timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu_pars(mu_parsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tp_pars(tp_parsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type m_pars(m_parsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ti_pars(ti_parsSEXP);
    Rcpp::traits::input_parameter< double >::type y0b(y0bSEXP);
    Rcpp::traits::input_parameter< double >::type lower_titre_bound(lower_titre_boundSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type times(timesSEXP);
    rcpp_result_gen = Rcpp::wrap(individual_sim(mu_pars, tp_pars, m_pars, ti_pars, y0b, lower_titre_bound, times));
    return rcpp_result_gen;
END_RCPP
}
// multiple_sim
NumericMatrix multiple_sim(NumericVector ti_pars, NumericVector y0s, NumericMatrix mu_pars, NumericMatrix tp_pars, NumericMatrix m_pars, NumericVector times);
RcppExport SEXP _serosim2_multiple_sim(SEXP ti_parsSEXP, SEXP y0sSEXP, SEXP mu_parsSEXP, SEXP tp_parsSEXP, SEXP m_parsSEXP, SEXP timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ti_pars(ti_parsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y0s(y0sSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mu_pars(mu_parsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type tp_pars(tp_parsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m_pars(m_parsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type times(timesSEXP);
    rcpp_result_gen = Rcpp::wrap(multiple_sim(ti_pars, y0s, mu_pars, tp_pars, m_pars, times));
    return rcpp_result_gen;
END_RCPP
}
// obs_error
double obs_error(int actual, int obs, double S, double EA);
RcppExport SEXP _serosim2_obs_error(SEXP actualSEXP, SEXP obsSEXP, SEXP SSEXP, SEXP EASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type actual(actualSEXP);
    Rcpp::traits::input_parameter< int >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< double >::type S(SSEXP);
    Rcpp::traits::input_parameter< double >::type EA(EASEXP);
    rcpp_result_gen = Rcpp::wrap(obs_error(actual, obs, S, EA));
    return rcpp_result_gen;
END_RCPP
}
// posterior_1
double posterior_1(NumericVector ti_pars, NumericVector y0s, double mu, NumericMatrix cr_pars, NumericMatrix tp_pars, NumericMatrix m_pars, NumericMatrix data, NumericVector error_pars, double mu_pop, double mu_pop_sigma);
RcppExport SEXP _serosim2_posterior_1(SEXP ti_parsSEXP, SEXP y0sSEXP, SEXP muSEXP, SEXP cr_parsSEXP, SEXP tp_parsSEXP, SEXP m_parsSEXP, SEXP dataSEXP, SEXP error_parsSEXP, SEXP mu_popSEXP, SEXP mu_pop_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ti_pars(ti_parsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y0s(y0sSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cr_pars(cr_parsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type tp_pars(tp_parsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m_pars(m_parsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type error_pars(error_parsSEXP);
    Rcpp::traits::input_parameter< double >::type mu_pop(mu_popSEXP);
    Rcpp::traits::input_parameter< double >::type mu_pop_sigma(mu_pop_sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(posterior_1(ti_pars, y0s, mu, cr_pars, tp_pars, m_pars, data, error_pars, mu_pop, mu_pop_sigma));
    return rcpp_result_gen;
END_RCPP
}
// posterior_2
double posterior_2(NumericVector ti_pars, NumericVector y0s, double mu, NumericMatrix cr_pars, NumericMatrix tp_pars, NumericMatrix m_pars, NumericMatrix data, double pop_sigma, double mu_pop, double mu_pop_sigma, NumericVector cr_pop, NumericVector cr_sigmas);
RcppExport SEXP _serosim2_posterior_2(SEXP ti_parsSEXP, SEXP y0sSEXP, SEXP muSEXP, SEXP cr_parsSEXP, SEXP tp_parsSEXP, SEXP m_parsSEXP, SEXP dataSEXP, SEXP pop_sigmaSEXP, SEXP mu_popSEXP, SEXP mu_pop_sigmaSEXP, SEXP cr_popSEXP, SEXP cr_sigmasSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ti_pars(ti_parsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y0s(y0sSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cr_pars(cr_parsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type tp_pars(tp_parsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m_pars(m_parsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type pop_sigma(pop_sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type mu_pop(mu_popSEXP);
    Rcpp::traits::input_parameter< double >::type mu_pop_sigma(mu_pop_sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cr_pop(cr_popSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cr_sigmas(cr_sigmasSEXP);
    rcpp_result_gen = Rcpp::wrap(posterior_2(ti_pars, y0s, mu, cr_pars, tp_pars, m_pars, data, pop_sigma, mu_pop, mu_pop_sigma, cr_pop, cr_sigmas));
    return rcpp_result_gen;
END_RCPP
}
// posterior_3
double posterior_3(NumericVector ti_pars, NumericVector y0s, double mu, NumericMatrix cr_pars, NumericMatrix tp_pars, NumericMatrix m_pars, NumericMatrix data, double pop_sigma, double mu_pop, double mu_pop_sigma, NumericVector cr_pop, NumericVector m_pop, NumericVector cr_sigmas, NumericVector m_sigmas);
RcppExport SEXP _serosim2_posterior_3(SEXP ti_parsSEXP, SEXP y0sSEXP, SEXP muSEXP, SEXP cr_parsSEXP, SEXP tp_parsSEXP, SEXP m_parsSEXP, SEXP dataSEXP, SEXP pop_sigmaSEXP, SEXP mu_popSEXP, SEXP mu_pop_sigmaSEXP, SEXP cr_popSEXP, SEXP m_popSEXP, SEXP cr_sigmasSEXP, SEXP m_sigmasSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ti_pars(ti_parsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y0s(y0sSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cr_pars(cr_parsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type tp_pars(tp_parsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m_pars(m_parsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type pop_sigma(pop_sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type mu_pop(mu_popSEXP);
    Rcpp::traits::input_parameter< double >::type mu_pop_sigma(mu_pop_sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cr_pop(cr_popSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type m_pop(m_popSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cr_sigmas(cr_sigmasSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type m_sigmas(m_sigmasSEXP);
    rcpp_result_gen = Rcpp::wrap(posterior_3(ti_pars, y0s, mu, cr_pars, tp_pars, m_pars, data, pop_sigma, mu_pop, mu_pop_sigma, cr_pop, m_pop, cr_sigmas, m_sigmas));
    return rcpp_result_gen;
END_RCPP
}
// posterior_4
double posterior_4(NumericVector ti_pars, NumericVector y0s, double mu, NumericMatrix cr_pars, NumericMatrix tp_pars, NumericMatrix m_pars, NumericMatrix data, NumericVector error_pars, double mu_pop, double mu_pop_sigma);
RcppExport SEXP _serosim2_posterior_4(SEXP ti_parsSEXP, SEXP y0sSEXP, SEXP muSEXP, SEXP cr_parsSEXP, SEXP tp_parsSEXP, SEXP m_parsSEXP, SEXP dataSEXP, SEXP error_parsSEXP, SEXP mu_popSEXP, SEXP mu_pop_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ti_pars(ti_parsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y0s(y0sSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cr_pars(cr_parsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type tp_pars(tp_parsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m_pars(m_parsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type error_pars(error_parsSEXP);
    Rcpp::traits::input_parameter< double >::type mu_pop(mu_popSEXP);
    Rcpp::traits::input_parameter< double >::type mu_pop_sigma(mu_pop_sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(posterior_4(ti_pars, y0s, mu, cr_pars, tp_pars, m_pars, data, error_pars, mu_pop, mu_pop_sigma));
    return rcpp_result_gen;
END_RCPP
}
// toUnitScale
double toUnitScale(double x, double min, double max);
RcppExport SEXP _serosim2_toUnitScale(SEXP xSEXP, SEXP minSEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type min(minSEXP);
    Rcpp::traits::input_parameter< double >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(toUnitScale(x, min, max));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_serosim2_likelihood_titre", (DL_FUNC) &_serosim2_likelihood_titre, 4},
    {"_serosim2_vector_sort", (DL_FUNC) &_serosim2_vector_sort, 1},
    {"_serosim2_vector_order", (DL_FUNC) &_serosim2_vector_order, 1},
    {"_serosim2_individual_sim", (DL_FUNC) &_serosim2_individual_sim, 7},
    {"_serosim2_multiple_sim", (DL_FUNC) &_serosim2_multiple_sim, 6},
    {"_serosim2_obs_error", (DL_FUNC) &_serosim2_obs_error, 4},
    {"_serosim2_posterior_1", (DL_FUNC) &_serosim2_posterior_1, 10},
    {"_serosim2_posterior_2", (DL_FUNC) &_serosim2_posterior_2, 12},
    {"_serosim2_posterior_3", (DL_FUNC) &_serosim2_posterior_3, 14},
    {"_serosim2_posterior_4", (DL_FUNC) &_serosim2_posterior_4, 10},
    {"_serosim2_toUnitScale", (DL_FUNC) &_serosim2_toUnitScale, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_serosim2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
