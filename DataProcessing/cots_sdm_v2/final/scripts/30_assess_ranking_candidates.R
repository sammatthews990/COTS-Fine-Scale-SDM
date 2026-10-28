# Question: do three fixed alternatives improve held-out search ranking?
# Method: forward GLMs with prior evidence, reef-held-out habitat means/tails;
# exact A-C site counts and matched retention, observed dives and bottom time.
# Validation: exact incumbent replay, frozen cohort and dates, paired reef-year
# bootstrap differences, recent-period check. No deployed product is replaced.
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))),'..','R','common.R'))
cfg <- initialise_final()
source(file.path(cfg$final_dir,'R','site_selection.R'))
source(file.path(cfg$final_dir,'R','ranking_candidates.R'))
cache <- file.path(cfg$dirs$cache,'ranking_candidates')
obj <- readRDS(file.path(cache,'features.rds'))
stopifnot(obj$signature$half_life_days==cfg$ranking_candidates$half_life_days)
ev <- obj$features;keys <- c('reef_name','year','site_row')
old <- readRDS(checked_input(cfg,'assessment_forward'))
cohort <- unique(old[c('reef_name','year')])
raw_ev <- readRDS(checked_input(cfg,'assessment_evidence')) |>
 dplyr::filter(.data$window=='last_season') |>
 dplyr::select(dplyr::all_of(keys),'site_problematic_prev','site_quiet_prev',
  'cull_rate','cull_effort','manta_strong','manta_effort','reef_cull_prev','reef_manta_prev') |>
 dplyr::inner_join(readRDS(checked_input(cfg,'assessment_manta_sites')) |>
  dplyr::select(dplyr::all_of(keys),'rank_cull','rank_within','rank_ensemble'),by=keys)
coefs <- read.csv(checked_input(cfg,'assessment_forward_coefficients'))
habitats <- list(cull=list(h='rank_cull'),ensemble=list(h='rank_ensemble'),both=list(h='rank_cull',h2='rank_within'))
set.seed(cfg$hmsc$seed+19L)
replayed <- list();raw_logits <- list()
for(T in cfg$ranking_candidates$evaluation_years) {
 f <- raw_ev[raw_ev$year==T,]
 scores <- data.frame(f[,keys],habitat_only=f$rank_cull)
 for(hn in names(habitats)) for(wn in c('all','last5','last3')) {
  h <- habitats[[hn]]
  x <- model.matrix(delete.response(terms(site_evidence_formula(habitat2=!is.null(h$h2)))),site_evidence_terms(f,h$h,h$h2))
  b <- coefs[coefs$test_year==T & coefs$habitat==hn & coefs$window==wn,]
  scores[[paste(hn,wn,sep='_')]] <- as.vector(x[,b$term,drop=FALSE]%*%b$estimate)
 }
 raw_logits[[as.character(T)]] <- data.frame(f[,keys],p_archived=plogis(scores$cull_all))
 sn <- setdiff(names(scores),keys)
 replayed[[as.character(T)]] <- scores |> dplyr::group_by(.data$reef_name,.data$year) |>
  dplyr::mutate(dplyr::across(dplyr::all_of(sn),pct_rank)) |> dplyr::ungroup()
}
replayed <- dplyr::bind_rows(replayed)
check <- dplyr::inner_join(old,replayed,by=keys,suffix=c('_old','_new'),relationship='many-to-one')
stopifnot(nrow(check)==nrow(old))
for(s in sn) stopifnot(identical(check[[paste0(s,'_old')]],check[[paste0(s,'_new')]]))
log_message('Every archived forward rank reproduced exactly.')
full <- replayed |> dplyr::semi_join(cohort,by=c('reef_name','year')) |>
 dplyr::select(dplyr::all_of(keys),archived_rank='cull_all') |>
 dplyr::group_by(.data$reef_name,.data$year) |>
 dplyr::mutate(N=dplyr::n(),pos_archived=ifelse(.data$N>1,round(1+(1-.data$archived_rank)*(.data$N-1)),1)) |>
 dplyr::ungroup() |>
 dplyr::left_join(dplyr::bind_rows(raw_logits),by=keys)
training <- obj$dives |> dplyr::filter(!is.na(.data$site_row)) |>
 dplyr::inner_join(ev,by=keys,relationship='many-to-one',suffix=c('','_feature'))
fits <- coeff <- list()
set.seed(cfg$hmsc$seed+30L)
ev$tie <- runif(nrow(ev))
for(T in cfg$ranking_candidates$evaluation_years) {
 tr <- training[training$year<T,]
 test <- ev[ev$year==T,]
 scores <- test[,keys]
 for(model in candidate_model_names) {
  x <- candidate_terms(tr,model);x$problematic <- as.numeric(tr$outcome=='problematic')
  fit <- stats::glm(candidate_formula(model),family=binomial(),data=x)
  stopifnot(fit$converged)
  test_score <- as.numeric(predict(fit,candidate_terms(test,model),type='link'))
  stopifnot(all(is.finite(test_score)))
  scores[[paste0('p_',model)]] <- plogis(test_score)
  temp <- data.frame(test[,c('reef_name','year')],s=test_score,tie=test$tie) |>
   dplyr::group_by(.data$reef_name,.data$year) |>
   dplyr::mutate(pos=candidate_position(.data$s,.data$tie)) |> dplyr::ungroup()
  scores[[paste0('pos_',model)]] <- temp$pos
  coeff[[length(coeff)+1L]] <- data.frame(test_year=T,model,term=names(coef(fit)),estimate=unname(coef(fit)),
   training_dives=nrow(tr),training_max_year=max(tr$year),converged=fit$converged)
  fits[[paste(T,model)]] <- fit
 }
 full <- dplyr::left_join(full,scores,by=keys,suffix=c('','_duplicate'))
 # Merge yearly columns rather than retaining one mostly-NA copy per year.
 for(model in candidate_model_names) for(prefix in c('p_','pos_')) {
  n <- paste0(prefix,model);dup <- paste0(n,'_duplicate')
  if(dup %in% names(full)) {full[[n]] <- dplyr::coalesce(full[[n]],full[[dup]]);full[[dup]] <- NULL}
 }
 log_message('Forward year ',T,': ',nrow(tr),' training dives; ',nrow(test),' candidates scored.')
}
write_table(dplyr::bind_rows(coeff),cfg,'ranking_candidate_forward_coefficients.csv')
full <- full |> dplyr::left_join(ev[,c(keys,'h_mean','h_q50','h_q75','h_q90','h_q95')],by=keys)
d <- dplyr::inner_join(obj$dives,full,by=keys,relationship='many-to-one')
stopifnot(nrow(d)==nrow(old),!anyDuplicated(full[keys]))
d$band <- ifelse(d$cpue>.04,'high',ifelse(d$cpue>=.02,'intermediate',ifelse(d$total_cots>0,'low','ghost')))
full$group <- paste(full$reef_name,full$year,sep='||');d$group <- paste(d$reef_name,d$year,sep='||')
groups <- sort(unique(full$group));ng <- length(groups)
full$g <- match(full$group,groups);d$g <- match(d$group,groups)
group_year <- as.integer(sub('.*\\|\\|','',groups))
ns <- tabulate(full$g,nbins=ng)
sum_groups <- function(x) {out <- rowsum(x*1,factor(d$g,levels=seq_len(ng)),reorder=TRUE);stopifnot(nrow(out)==ng);out}
fractions <- seq(0,1,.01);a80 <- which(abs(fractions-.8)<1e-10)
den <- lapply(c('high','intermediate','low','ghost'),function(b) as.numeric(sum_groups(as.numeric(d$band==b))))
names(den) <- c('high','intermediate','low','ghost')
den$dives <- as.numeric(sum_groups(rep(1,nrow(d))));den$minutes <- as.numeric(sum_groups(d$bottom_time))
models <- c('archived',candidate_model_names,'random')
counts <- auc_parts <- list()
auc_piece <- function(p,y,include=rep(TRUE,nrow(d))) {
 out <- matrix(0,ng,2)
 for(g in seq_len(ng)) {
  ii <- which(d$g==g & include);yy <- y[ii];n1 <- sum(yy);n0 <- sum(!yy)
  if(n1>0 && n0>0) out[g,] <- c(sum(rank(p[ii])[yy])-n1*(n1+1)/2,n1*n0)
 }
 out
}
for(model in models) {
 keep <- if(model=='random') sapply(fractions,function(a) floor(a*d$N+1e-9)/d$N) else
  sapply(fractions,function(a) d[[paste0('pos_',model)]]<=floor(a*d$N+1e-9))
 z <- lapply(c('high','intermediate','low','ghost'),function(b) sum_groups(keep*(d$band==b)))
 names(z) <- c('high','intermediate','low','ghost')
 z$dives <- sum_groups(keep);z$minutes <- sum_groups(keep*d$bottom_time)
 z$sites <- sapply(fractions,function(a) floor(a*ns+1e-9))
 counts[[model]] <- z
 if(model!='random') {
  p <- d[[paste0('p_',model)]];p <- pmin(1-1e-8,pmax(1e-8,p));y <- d$cpue>=.02
  auc_parts[[model]] <- list(high_all=auc_piece(p,d$cpue>.04),
   high_quiet=auc_piece(p,d$cpue>.04,d$cpue>.04 | d$cpue<.02),
   useful=auc_piece(p,y),brier=as.numeric(sum_groups((p-y)^2)),
   logloss=as.numeric(sum_groups(-y*log(p)-(1-y)*log(1-p))))
 }
}
# Habitat-only scores isolate whether aggregation obscured useful habitat.
for(h in c('h_mean','h_q50','h_q75','h_q90','h_q95')) {
 auc_parts[[h]] <- list(high_all=auc_piece(d[[h]],d$cpue>.04),
  useful=auc_piece(d[[h]],d$cpue>=.02))
}
periods <- list(`2019-2026`=2019:2026,`2024-2026`=2024:2026)
periods <- c(periods,setNames(as.list(2019:2026),as.character(2019:2026)))
all_values <- rows <- curves <- list()
interval <- function(v) c(value=v[1],low=unname(quantile(v[-1],.025,na.rm=TRUE)),high=unname(quantile(v[-1],.975,na.rm=TRUE)))
add <- function(model,period,rule,target,metric,v) {
 key <- paste(model,period,rule,target,metric,sep='::')
 all_values[[key]] <<- v
 rows[[length(rows)+1L]] <<- data.frame(model,period,rule,target,metric,t(interval(v)),check.names=FALSE)
}
for(period in names(periods)) {
 ii <- which(group_year %in% periods[[period]]);n <- length(ii)
 set.seed(cfg$hmsc$seed+30L+min(periods[[period]]))
 W <- cbind(rep(1,n),replicate(cfg$ranking_candidates$bootstrap_reps,tabulate(sample.int(n,n,replace=TRUE),nbins=n)))
 D <- lapply(den,function(x) as.vector(crossprod(W,x[ii])))
 site_den <- as.vector(crossprod(W,ns[ii]))
 ratio <- function(x,total) crossprod(W,x[ii,,drop=FALSE])/total
 out <- list()
 for(model in models) {
  cc <- counts[[model]]
  v <- list(sites_selected=ratio(cc$sites,site_den),high_retained=ratio(cc$high,D$high),
   intermediate_retained=ratio(cc$intermediate,D$intermediate),
   useful_retained=ratio(cc$high+cc$intermediate,D$high+D$intermediate),
   low_avoided=1-ratio(cc$low,D$low),ghost_avoided=1-ratio(cc$ghost,D$ghost),
   dives_avoided=1-ratio(cc$dives,D$dives),bottom_time_avoided=1-ratio(cc$minutes,D$minutes))
  out[[model]] <- v
  for(metric in names(v)) {
   add(model,period,'A-C',.8,metric,v[[metric]][,a80])
   for(target in c(.8,.9,.95,.99)) add(model,period,'matched_high',target,metric,
    first_crossing(v$high_retained,v[[metric]],target))
  }
  quiet <- (D$low+D$ghost)-crossprod(W,(cc$low+cc$ghost)[ii,a80])
  deferred <- D$dives-crossprod(W,cc$dives[ii,a80])
  add(model,period,'A-C',.8,'quiet_share_deferred',as.vector(quiet/deferred))
  top20 <- which(abs(fractions-.2)<1e-10)
  random_hits <- crossprod(W,counts$random$high[ii,top20])
  add(model,period,'discrimination',.2,'top20_lift',as.vector(crossprod(W,cc$high[ii,top20])/random_hits))
  if(model!='random') {
   ap <- auc_parts[[model]]
   for(name in c('high_all','high_quiet','useful')) {
    a <- ap[[name]]
    add(model,period,'discrimination',0,paste0('auc_',name),as.vector(crossprod(W,a[ii,1])/crossprod(W,a[ii,2])))
   }
   for(name in c('brier','logloss')) add(model,period,'discrimination',0,name,
    as.vector(crossprod(W,ap[[name]][ii]))/D$dives)
  }
  curves[[length(curves)+1L]] <- data.frame(model,period,site_fraction=fractions,
   sites_selected=v$sites_selected[1,],high_retained=v$high_retained[1,],
   intermediate_retained=v$intermediate_retained[1,],dives_avoided=v$dives_avoided[1,],
   bottom_time_avoided=v$bottom_time_avoided[1,])
 }
 for(model in models) for(budget in c('dives','bottom_time')) {
  target <- 1-out$archived[[paste0(budget,'_avoided')]][,a80]
  x <- 1-out[[model]][[paste0(budget,'_avoided')]]
  for(metric in names(out[[model]])) add(model,period,paste0('matched_',budget),.8,metric,
   first_crossing(x,out[[model]][[metric]],target))
 }
 for(h in c('h_mean','h_q50','h_q75','h_q90','h_q95')) for(name in c('high_all','useful')) {
  a <- auc_parts[[h]][[name]]
  add(h,period,'habitat_only',0,paste0('auc_',name),as.vector(crossprod(W,a[ii,1])/crossprod(W,a[ii,2])))
 }
}
results <- dplyr::bind_rows(rows)
write_table(results,cfg,'ranking_candidate_assessment.csv')
write_table(dplyr::bind_rows(curves),cfg,'ranking_candidate_curves.csv')
differences <- list()
for(i in seq_len(nrow(results))) {
 r <- results[i,];if(!r$model %in% c(candidate_model_names,'h_q50','h_q75','h_q90','h_q95')) next
 references <- if(grepl('^h_',r$model)) 'h_mean' else c('archived','mean_control','random')
 for(reference in references) {
  if(reference==r$model) next
  key <- paste(r$model,r$period,r$rule,r$target,r$metric,sep='::')
  ref <- paste(reference,r$period,r$rule,r$target,r$metric,sep='::')
  if(!ref %in% names(all_values)) next
  delta <- all_values[[key]]-all_values[[ref]]
  differences[[length(differences)+1L]] <- data.frame(r[,c('model','period','rule','target','metric')],reference,t(interval(delta)),check.names=FALSE)
 }
}
differences <- dplyr::bind_rows(differences)
write_table(differences,cfg,'ranking_candidate_differences.csv')
policy <- read.csv(file.path(cfg$dirs$tables,'operational_priority_cutoff.csv'))
for(metric in c('sites_selected','high_retained','low_avoided','ghost_avoided')) {
 a <- results$value[results$model=='archived' & results$period=='2019-2026' & results$rule=='A-C' & results$metric==metric]
 b <- policy$value[policy$score=='combined' & policy$metric==metric]
 stopifnot(abs(a-b)<1e-12)
}
saveRDS(list(full=full,dives=d,results=results,differences=differences,counts=counts,denominators=den,
 groups=groups,group_year=group_year,curves=dplyr::bind_rows(curves),fits=fits,bootstrap_values=all_values),
 file.path(cache,'assessment.rds'))
show <- results[results$period=='2019-2026' & results$rule=='matched_high' & results$target==.9 &
 results$metric %in% c('sites_selected','dives_avoided','bottom_time_avoided','intermediate_retained'),]
print(show,row.names=FALSE)
log_message('Three candidate assessments complete. Incumbent product is unchanged.')
