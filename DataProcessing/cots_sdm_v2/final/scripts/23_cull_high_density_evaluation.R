# Do existing held-out site rankings retain dives above 0.04 CPUE while avoiding
# empty and low-density dives? Re-score fixed rankings, without refitting or
# changing the 2027 product. Compare at identical site effort and at matched
# high-density retention; use paired reef-year bootstrap intervals (500 draws).
# Spatial rankings are from 15-17; forward evidence rankings are from 19.
# Habitat fits in that forward comparison are spatially, not temporally, held out.
source(file.path(dirname(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE))), '..', 'R', 'common.R'))
cfg <- initialise_final()
source(file.path(cfg$final_dir,'R/cull_sites.R'))
set.seed(cfg$hmsc$seed + 23L)
keys <- c('reef_name','year','site_row')
pin_names <- grep('^assessment_', names(cfg$inputs), value=TRUE)
paths <- stats::setNames(lapply(pin_names, function(n) checked_input(cfg,n)), pin_names)
write_table(data.frame(input=pin_names, sha256=vapply(paths,sha256_file,character(1))),cfg,'cull_high_density_input_manifest.csv')
# The matching helper reads the same checked files, in their original order.
stopifnot(normalizePath(paths$assessment_dives)==normalizePath(file.path(cfg$exploration_dir,'artifacts/cull_observations_clean.rds')))
sites <- read_designated_sites(cfg)
dives <- matched_cull_dives(cfg,sites)
stopifnot(all(dives$total_cots>=0), all(is.finite(dives$cpue)))
write_table(data.frame(valid_dives=nrow(dives), all_record_reefs=dplyr::n_distinct(dives$reef_name),
  matched_reefs=dplyr::n_distinct(dives$reef_name[!is.na(dives$site_row)]),
  matched_sites=dplyr::n_distinct(na.omit(dives$site_row)), first_year=min(dives$year),last_year=max(dives$year)),
  cfg,'cull_record_scope.csv')
dives$band <- ifelse(dives$total_cots==0,'ghost',ifelse(dives$cpue<0.02,'low',ifelse(dives$cpue<=0.04,'intermediate','high')))
stopifnot(all(dives$band[dives$cpue==0.02]=='intermediate'),all(dives$band[dives$cpue==0.04]=='intermediate'))
e <- readRDS(paths$assessment_evidence)
e <- e[e$window=='last_season',]
m <- readRDS(paths$assessment_manta_sites)
sp <- e[,c(keys,'rank_cv','rank_evidence','rank_combined')]
names(sp)[4:6] <- c('habitat','evidence','combined')
sp <- dplyr::inner_join(sp,m[,c(keys,'rank_within','rank_ensemble')],by=keys,relationship='one-to-one')
names(sp)[names(sp)=='rank_within'] <- 'manta_habitat'
names(sp)[names(sp)=='rank_ensemble'] <- 'habitat_ensemble'
fr <- readRDS(paths$assessment_forward)
score_cols <- c('habitat_only','cull_all','ensemble_all','both_all')
fw <- unique(fr[,c(keys,score_cols)])
stopifnot(!anyDuplicated(fw[keys]))
names(fw)[4:7] <- c('habitat','combined','ensemble_combined','two_habitat_combined')
matched <- dives[!is.na(dives$site_row),]
fd <- dplyr::inner_join(matched,fw,by=keys,relationship='many-to-one')
# Scores were saved once per dive by 19; prove that the join recovers that exact
# evaluation sample, avoiding a many-to-many join of repeated dive-site keys.
check <- dplyr::full_join(dplyr::count(fr,reef_name,year,site_row,outcome,name='old_n'),
 dplyr::count(fd,reef_name,year,site_row,outcome,name='new_n'),by=c(keys,'outcome'))
stopifnot(nrow(fd)==nrow(fr),all(check$old_n==check$new_n),!anyNA(check))

fractions <- round(seq(0.05,1,0.05),2)
targets <- c(0.8,0.9,0.95,0.99)
bands <- c('ghost','low','intermediate','high')
reps <- 500L
fixed_rows <- matched_rows <- difference_rows <- composition_rows <- list()

run_set <- function(candidates,design,period,years,effort_candidates=candidates) {
 score_names <- setdiff(names(candidates),keys)
 stopifnot(all(is.finite(as.matrix(candidates[,score_names]))))
 candidates <- candidates[candidates$year %in% years,]
 d <- dplyr::inner_join(matched,candidates,by=keys,relationship='many-to-one')
 d$reef_year <- paste(d$reef_name,d$year,sep='||')
 candidates <- effort_candidates[effort_candidates$year %in% years,]
 candidates$reef_year <- paste(candidates$reef_name,candidates$year,sep='||')
 candidates <- candidates[candidates$reef_year %in% d$reef_year,]
 groups <- sort(unique(d$reef_year)); ng <- length(groups)
 for (band in bands) {
  x <- d[d$band==band,]
  composition_rows[[length(composition_rows)+1L]] <<- data.frame(design,period,band,dives=nrow(x),
    reefs=dplyr::n_distinct(x$reef_name),reef_years=dplyr::n_distinct(x$reef_year))
 }
 totals <- sapply(bands,function(b) as.numeric(table(factor(d$reef_year[d$band==b],levels=groups))))
 ns <- as.numeric(table(factor(candidates$reef_year,levels=groups)))
 for(sc in score_names) {
  positions <- d[[sc]]*pmax(ns[match(d$reef_year,groups)]-1,1)
  # Script 15 used mean ranks for ties; 16-19 used random tie breaking.
  if(design=='forward') stopifnot(all(abs(positions-round(positions))<1e-7))
  stopifnot(all(d[[sc]]>=0 & d[[sc]]<=1))
 }
 # Every draw uses the same reef-year multiplicities for every ranking/outcome.
 weights <- cbind(rep(1,ng),replicate(reps,tabulate(sample.int(ng,ng,replace=TRUE),nbins=ng)))
 den <- t(weights)%*%totals
 nsite <- as.vector(t(weights)%*%ns)
 curves <- list()
 for (sc in c(score_names,'random')) {
  ret <- avoid_low <- avoid_ghost <- intermediate <- searched <- matrix(NA_real_,reps+1L,length(fractions)+1L)
  ret[,1] <- intermediate[,1] <- searched[,1] <- 0
  avoid_low[,1] <- avoid_ghost[,1] <- 1
  for (j in seq_along(fractions)) {
   a <- fractions[j]
   # Random expectation uses the exact count selected on each reef by the
   # reference ranking, so its effort matches even on reefs with few sites.
   # Count spatial selections directly, including tied habitat ranks. Forward
   # ranks span all designated candidates with random tie breaking, so their
   # exact counts follow from N even though only dived-site scores were saved.
   if(design=='spatial') {
    ref <- if(sc=='random') 'combined' else sc
    keep_s <- candidates[[ref]]>=1-a-1e-9
    k <- as.numeric(table(factor(candidates$reef_year[keep_s],levels=groups)))
   } else k <- vapply(ns,function(n) sum(seq(0,1,length.out=n)>=1-a-1e-9),numeric(1))
   if(sc=='random') counts <- totals*(k/ns) else {
    keep <- d[[sc]]>=1-a-1e-9
    counts <- sapply(bands,function(b) as.numeric(table(factor(d$reef_year[keep & d$band==b],levels=groups))))
   }
   num <- t(weights)%*%counts
   ret[,j+1] <- num[,4]/den[,4]
   avoid_low[,j+1] <- 1-num[,2]/den[,2]
   avoid_ghost[,j+1] <- 1-num[,1]/den[,1]
   intermediate[,j+1] <- num[,3]/den[,3]
   searched[,j+1] <- as.vector(t(weights)%*%k)/nsite
  }
  stopifnot(all(is.finite(ret)),all(abs(ret[,ncol(ret)]-1)<1e-10),
    all(abs(searched[,ncol(searched)]-1)<1e-10),all(abs(avoid_low[,ncol(avoid_low)])<1e-10),
    all(abs(avoid_ghost[,ncol(avoid_ghost)])<1e-10),all(apply(ret,1,function(x) all(diff(x)>=-1e-10))))
  curves[[sc]] <- list(high_retained=ret,high_missed=1-ret,low_avoided=avoid_low,ghost_avoided=avoid_ghost,
    intermediate_retained=intermediate,sites_searched=searched)
 }
 interval <- function(x) c(value=x[1],low=unname(quantile(x[-1],.025,na.rm=TRUE)),high=unname(quantile(x[-1],.975,na.rm=TRUE)))
 for(sc in names(curves)) for(metric in names(curves[[sc]])) for(j in seq_along(fractions)) {
  v <- interval(curves[[sc]][[metric]][,j+1])
  fixed_rows[[length(fixed_rows)+1L]] <<- data.frame(design,period,score=sc,nominal_fraction=fractions[j],metric,
    value=v[1],low=v[2],high=v[3],dives=nrow(d),reef_years=ng)
 }
 # First crossing gives the minimum site effort attaining the target. Linear
 # interpolation only between neighbouring effort levels; never average plateaus.
 matched_curves <- list()
 for(sc in names(curves)) for(target in targets) {
  out <- sapply(names(curves[[sc]]),function(metric) vapply(seq_len(reps+1),function(b) {
   x <- curves[[sc]]$high_retained[b,]; y <- curves[[sc]][[metric]][b,]
   if(any(!is.finite(x))) return(NA_real_)
   j <- which(x>=target)[1]
   if(x[j]==target) y[j] else y[j-1]+(y[j]-y[j-1])*(target-x[j-1])/(x[j]-x[j-1])
  },numeric(1)))
  matched_curves[[paste(sc,target)]] <- out
  for(metric in colnames(out)) {
   v <- interval(out[,metric])
   matched_rows[[length(matched_rows)+1L]] <<- data.frame(design,period,score=sc,target,metric,value=v[1],low=v[2],high=v[3],dives=nrow(d),reef_years=ng)
  }
 }
 for(target in targets) for(baseline in c('habitat','random')) for(metric in c('sites_searched','low_avoided','ghost_avoided','intermediate_retained')) {
  delta <- matched_curves[[paste('combined',target)]][,metric]-matched_curves[[paste(baseline,target)]][,metric]
  v <- interval(delta)
  difference_rows[[length(difference_rows)+1L]] <<- data.frame(design,period,score='combined',baseline,target,metric,value=v[1],low=v[2],high=v[3])
 }
}
run_set(sp,'spatial','2018-2026',2018:2026)
run_set(sp,'spatial','2022-2026',2022:2026)
run_set(fw,'forward','2019-2026',2019:2026,sp)
run_set(fw,'forward','2022-2026',2022:2026,sp)
for(y in 2019:2026) run_set(fw,'forward',as.character(y),y,sp)
fixed <- dplyr::bind_rows(fixed_rows); mt <- dplyr::bind_rows(matched_rows)
stopifnot(all(is.finite(fixed$value)),all(is.finite(mt$value)),all(fixed$value>=-1e-10 & fixed$value<=1+1e-10))
write_table(fixed,cfg,'cull_high_density_fixed.csv')
write_table(mt,cfg,'cull_high_density_matched.csv')
write_table(dplyr::bind_rows(difference_rows),cfg,'cull_high_density_differences.csv')
write_table(dplyr::bind_rows(composition_rows),cfg,'cull_high_density_composition.csv')
write_table(fixed[fixed$design=='forward' & fixed$period=='2019-2026' & fixed$score=='combined' &
  fixed$nominal_fraction %in% c(.2,.5,.8,1),],cfg,'cull_high_density_priority_classes.csv')
log_message('Completed fixed-ranking assessment; no refits or deployed-priority changes.')
