# Helpers for strictly forward (train < T, test = T) validation of habitat,
# reef-context and history-kernel scores.

context_term_names <- c(
  "edge_signed_log", "windward_position", "seaward_position", "log_upwind_km",
  "log_gap_km", "log_reef_area_km2", "log_current_speed", "flow_steadiness"
)

# Transform raw reef-context values into model terms. `impute` carries the
# training medians so prediction uses the same fills.
context_terms <- function(raw, impute = NULL) {
  terms <- data.frame(
    edge_signed_log = sign(raw$dist_to_reef_edge_m) * log1p(abs(raw$dist_to_reef_edge_m) / 30),
    windward_position = raw$windward_position,
    seaward_position = raw$seaward_position,
    log_upwind_km = log(raw$upwind_open_water_km),
    log_gap_km = log1p(raw$nearest_other_reef_km),
    log_reef_area_km2 = raw$log_reef_area_km2,
    log_current_speed = log(pmax(raw$current_speed, 1e-3)),
    flow_steadiness = raw$flow_steadiness
  )
  if (is.null(impute)) impute <- vapply(terms, stats::median, numeric(1), na.rm = TRUE)
  for (j in names(terms)) terms[[j]][!is.finite(terms[[j]])] <- impute[[j]]
  attr(terms, "impute") <- impute
  terms
}

# Kernel-weighted prior high-density rate around each query point.
#   weight = exp(-d^2 / (2 sigma^2)) * exp(-(age) / tau), truncated at 3 sigma.
# Returns rate H and total weight W (W = 0 means no usable history).
history_kernel <- function(prior_xy, prior_y, prior_age, query_xy, sigma, tau, k = 400L) {
  n <- nrow(query_xy)
  if (!nrow(prior_xy)) return(data.frame(H = rep(NA_real_, n), W = rep(0, n)))
  nn <- FNN::get.knnx(prior_xy, query_xy, k = min(k, nrow(prior_xy)))
  w <- exp(-nn$nn.dist^2 / (2 * sigma^2)) * (nn$nn.dist <= 3 * sigma)
  if (is.finite(tau)) w <- w * exp(-matrix(prior_age[nn$nn.index], nrow = n) / tau)
  y <- matrix(prior_y[nn$nn.index], nrow = n)
  W <- rowSums(w)
  data.frame(H = ifelse(W > 0, rowSums(w * y) / W, NA_real_), W = W)
}

# Habitat prior updated by history: (kappa * p_habitat + sum(w y)) / (kappa + W).
combine_history <- function(p_habitat, H, W, kappa) {
  out <- (kappa * p_habitat + ifelse(W > 0, H * W, 0)) / (kappa + W)
  out[!is.finite(p_habitat)] <- H[!is.finite(p_habitat)]
  out
}
