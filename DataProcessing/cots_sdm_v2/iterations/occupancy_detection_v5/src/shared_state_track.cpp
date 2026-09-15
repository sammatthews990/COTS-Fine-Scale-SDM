#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_MATRIX(X);
  DATA_IVECTOR(has_point_observation);

  DATA_IVECTOR(direct_site);
  DATA_IVECTOR(direct_method);
  DATA_IVECTOR(direct_y);
  DATA_VECTOR(direct_log_effort);
  DATA_VECTOR(direct_post_cull);
  DATA_VECTOR(direct_log_days_after_cull);

  DATA_IVECTOR(scar_site);
  DATA_IVECTOR(scar_method);
  DATA_IVECTOR(scar_y);
  DATA_VECTOR(scar_prior_direct_known);
  DATA_VECTOR(scar_log_days_since_direct);
  DATA_VECTOR(scar_post_cull);

  DATA_IVECTOR(manta_group_start);
  DATA_IVECTOR(manta_group_type);
  DATA_IVECTOR(manta_link_site);
  DATA_VECTOR(manta_link_log_effort);
  DATA_VECTOR(manta_group_weight);

  PARAMETER_VECTOR(beta_state);
  PARAMETER_VECTOR(alpha_direct);
  PARAMETER_VECTOR(effort_direct);
  PARAMETER(depletion_post);
  PARAMETER(depletion_log_days);
  PARAMETER_VECTOR(alpha_scar);
  PARAMETER_VECTOR(effort_scar);
  PARAMETER_VECTOR(logit_fp_scar);
  PARAMETER(scar_prior_effect);
  PARAMETER(scar_decay);
  PARAMETER(scar_post_cull_effect);
  PARAMETER(alpha_manta_direct);
  PARAMETER(effort_manta_direct);
  PARAMETER(alpha_manta_scar);
  PARAMETER(effort_manta_scar);
  PARAMETER(logit_fp_manta_scar);

  const Type eps = Type(1e-10);
  int n_site = X.rows();
  vector<Type> psi = (X * beta_state).unaryExpr([](Type x) {
    return invlogit(x);
  });
  vector<Type> log_lik_z1(n_site);
  vector<Type> log_lik_z0(n_site);
  log_lik_z1.setZero();
  log_lik_z0.setZero();
  Type nll = 0;

  for (int j = 0; j < direct_y.size(); j++) {
    int site = direct_site(j);
    int method = direct_method(j);
    Type eta = alpha_direct(method) + effort_direct(method) * direct_log_effort(j);
    if (method == 0) {
      eta += depletion_post * direct_post_cull(j) +
        depletion_log_days * direct_log_days_after_cull(j);
    }
    Type p1 = invlogit(eta);
    Type p0 = eps;
    log_lik_z1(site) += dbinom(Type(direct_y(j)), Type(1), p1, true);
    log_lik_z0(site) += dbinom(Type(direct_y(j)), Type(1), p0, true);
  }

  for (int j = 0; j < scar_y.size(); j++) {
    int site = scar_site(j);
    int method = scar_method(j);
    Type eta = alpha_scar(method) + effort_scar(method) * direct_log_effort(0) * Type(0);
    eta += scar_prior_effect * scar_prior_direct_known(j) +
      scar_decay * scar_log_days_since_direct(j) +
      scar_post_cull_effect * scar_post_cull(j);
    Type p1 = invlogit(eta);
    Type p0 = invlogit(logit_fp_scar(method));
    log_lik_z1(site) += dbinom(Type(scar_y(j)), Type(1), p1, true);
    log_lik_z0(site) += dbinom(Type(scar_y(j)), Type(1), p0, true);
  }

  for (int i = 0; i < n_site; i++) {
    if (has_point_observation(i) == 1) {
      Type occupied = log(psi(i) + eps) + log_lik_z1(i);
      Type unoccupied = log(Type(1) - psi(i) + eps) + log_lik_z0(i);
      nll -= logspace_add(occupied, unoccupied);
    }
  }

  for (int g = 0; g < manta_group_type.size(); g++) {
    Type log_no_detection = 0;
    for (int k = manta_group_start(g); k < manta_group_start(g + 1); k++) {
      int site = manta_link_site(k);
      Type p_direct = invlogit(alpha_manta_direct +
        effort_manta_direct * manta_link_log_effort(k));
      Type p_scar = invlogit(alpha_manta_scar +
        effort_manta_scar * manta_link_log_effort(k));
      Type fp_scar = invlogit(logit_fp_manta_scar);
      Type q_direct = psi(site) * p_direct + (Type(1) - psi(site)) * eps;
      Type q_scar = psi(site) * p_scar + (Type(1) - psi(site)) * fp_scar;
      Type q = q_direct;
      if (manta_group_type(g) == 2) q = q_scar;
      if (manta_group_type(g) == 0) {
        q = Type(1) - (Type(1) - q_direct) * (Type(1) - q_scar);
      }
      log_no_detection += log(Type(1) - q + eps);
    }
    Type group_log_lik = log_no_detection;
    if (manta_group_type(g) > 0) {
      group_log_lik = log(Type(1) - exp(log_no_detection) + eps);
    }
    nll -= manta_group_weight(g) * group_log_lik;
  }

  for (int k = 0; k < beta_state.size(); k++) nll -= dnorm(beta_state(k), Type(0), Type(2), true);
  for (int k = 0; k < alpha_direct.size(); k++) {
    nll -= dnorm(alpha_direct(k), Type(0), Type(2.5), true);
    nll -= dnorm(effort_direct(k), Type(0), Type(1), true);
  }
  for (int k = 0; k < alpha_scar.size(); k++) {
    nll -= dnorm(alpha_scar(k), Type(0), Type(2.5), true);
    nll -= dnorm(effort_scar(k), Type(0), Type(1), true);
    nll -= dnorm(logit_fp_scar(k), Type(-4), Type(1.5), true);
  }
  nll -= dnorm(depletion_post, Type(0), Type(1.5), true);
  nll -= dnorm(depletion_log_days, Type(0), Type(1), true);
  nll -= dnorm(scar_prior_effect, Type(0), Type(1.5), true);
  nll -= dnorm(scar_decay, Type(0), Type(1), true);
  nll -= dnorm(scar_post_cull_effect, Type(0), Type(1.5), true);
  nll -= dnorm(alpha_manta_direct, Type(0), Type(2.5), true);
  nll -= dnorm(effort_manta_direct, Type(0), Type(1), true);
  nll -= dnorm(alpha_manta_scar, Type(0), Type(2.5), true);
  nll -= dnorm(effort_manta_scar, Type(0), Type(1), true);
  nll -= dnorm(logit_fp_manta_scar, Type(-4), Type(1.5), true);

  ADREPORT(psi);
  REPORT(log_lik_z1);
  REPORT(log_lik_z0);
  return nll;
}
