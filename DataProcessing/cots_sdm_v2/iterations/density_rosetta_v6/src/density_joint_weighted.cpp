#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_MATRIX(X_cull);
  DATA_IVECTOR(cull_y);
  DATA_VECTOR(cull_effort);
  DATA_MATRIX(X_salad_link);
  DATA_IVECTOR(salad_link_track);
  DATA_VECTOR(salad_link_area_ha);
  DATA_IVECTOR(salad_y);
  DATA_SCALAR(calibration_alpha);
  DATA_SCALAR(calibration_beta);
  DATA_SCALAR(calibration_log_density_offset);
  DATA_SCALAR(cull_likelihood_weight);
  DATA_SCALAR(salad_likelihood_weight);

  PARAMETER_VECTOR(beta);
  PARAMETER(log_phi_cull);
  PARAMETER(log_phi_salad);

  Type nll = 0;
  Type phi_cull = exp(log_phi_cull);
  Type phi_salad = exp(log_phi_salad);

  if (cull_y.size() > 0) {
    vector<Type> eta_cull = X_cull * beta;
    for (int i = 0; i < cull_y.size(); i++) {
      Type log_rate = calibration_alpha + calibration_beta *
        (eta_cull(i) + calibration_log_density_offset);
      Type mu = cull_effort(i) * exp(log_rate);
      Type variance = mu + mu * mu / phi_cull;
      nll -= cull_likelihood_weight *
        dnbinom2(Type(cull_y(i)), mu, variance, true);
    }
  }

  if (salad_y.size() > 0) {
    vector<Type> eta_link = X_salad_link * beta;
    vector<Type> mu_salad(salad_y.size());
    mu_salad.setZero();
    for (int j = 0; j < salad_link_track.size(); j++) {
      mu_salad(salad_link_track(j)) += salad_link_area_ha(j) * exp(eta_link(j));
    }
    for (int i = 0; i < salad_y.size(); i++) {
      Type variance = mu_salad(i) + mu_salad(i) * mu_salad(i) / phi_salad;
      nll -= salad_likelihood_weight *
        dnbinom2(Type(salad_y(i)), mu_salad(i), variance, true);
    }
  }

  for (int k = 0; k < beta.size(); k++) {
    Type prior_sd = (k == 0) ? Type(5) : Type(1.5);
    nll -= dnorm(beta(k), Type(0), prior_sd, true);
  }
  nll -= dnorm(log_phi_cull, log(Type(5)), Type(1.5), true);
  nll -= dnorm(log_phi_salad, log(Type(5)), Type(1.5), true);

  ADREPORT(beta);
  REPORT(phi_cull);
  REPORT(phi_salad);
  return nll;
}
