# Candidate definitions fixed in report/ranking_alternatives_analysis_plan.md.
candidate_model_names <- c('mean_control','high_indicator','recent_cpue','hotspot_cpue')
candidate_terms <- function(f,model) {
 h <- if(model=='hotspot_cpue') 'h_q90' else 'h_mean'
 x <- site_evidence_terms(f,h)
 if(model!='mean_control') x$site_problematic_prev <- f$site_high_prev
 if(model %in% c('recent_cpue','hotspot_cpue')) {
  x$log_recent_cpue <- log1p(f$recent_cpue/.02)
  x$log_recent_weight <- log1p(f$recent_weight)
 }
 x
}
candidate_formula <- function(model) {
 terms <- c('habitat','site_problematic_prev','site_quiet_prev','cull_rate',
  'log_cull_effort','log_manta_strong','log_manta_effort')
 if(model %in% c('recent_cpue','hotspot_cpue')) terms <- c(terms,'log_recent_cpue','log_recent_weight')
 stats::reformulate(terms,response='problematic')
}
candidate_position <- function(score,tie) {
 out <- integer(length(score));out[order(-score,tie)] <- seq_along(score);out
}
first_crossing <- function(x,y,target) {
 vapply(seq_len(nrow(x)),function(i) {
  t <- if(length(target)==1) target else target[i]
  j <- which(x[i,]>=t-1e-12)[1]
  if(is.na(j)) return(NA_real_)
  if(j==1 || abs(x[i,j]-t)<1e-12) return(y[i,j])
  y[i,j-1]+(y[i,j]-y[i,j-1])*(t-x[i,j-1])/(x[i,j]-x[i,j-1])
 },numeric(1))
}
