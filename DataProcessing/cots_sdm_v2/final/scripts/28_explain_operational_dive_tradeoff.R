# Question: which observed dives would the exact A-C cutoff defer?
# Method: use frozen forward ranks and four exhaustive CPUE classes; compare
# counts and deferred-dive composition with exact same-site-count randomness.
# Validation: reproduce script 27's three published dive-outcome point estimates.
# Paired uncertainty resamples whole reef-years 500 times, matching script 27.
# No refitting or changes to the deployed product.
source('DataProcessing/cots_sdm_v2/final/R/common.R')
cfg <- initialise_final()
source(file.path(cfg$final_dir, 'R/cull_sites.R'))
e <- readRDS(checked_input(cfg,'assessment_evidence'))
e <- e[e$window=='last_season',]
f <- readRDS(checked_input(cfg,'assessment_forward'))
n <- dplyr::count(e,reef_name,year,name='N')
n <- dplyr::semi_join(n,unique(f[c('reef_name','year')]),by=c('reef_name','year'))
d <- matched_cull_dives(cfg,read_designated_sites(cfg))
d <- dplyr::inner_join(d,unique(f[c('reef_name','year','site_row','cull_all')]),
 by=c('reef_name','year','site_row'),relationship='many-to-one')
d <- dplyr::left_join(d,n,by=c('reef_name','year'))
stopifnot(nrow(d)==nrow(f))
d$position <- ifelse(d$N>1,round(1+(1-d$cull_all)*(d$N-1)),1)
d$deferred <- d$position/d$N>.8
d$random_deferred <- 1-floor(.8*d$N+1e-9)/d$N
d$band <- ifelse(d$cpue>.04,'high',ifelse(d$total_cots==0,'ghost',
 ifelse(d$cpue<.02,'low','intermediate')))
out <- dplyr::summarise(dplyr::group_by(d,.data$band),
 observed_dives=dplyr::n(),deferred_dives=sum(.data$deferred),
 random_deferred_dives=sum(.data$random_deferred),.groups='drop')
out$share_all_dives <- out$observed_dives/sum(out$observed_dives)
out$avoided <- out$deferred_dives/out$observed_dives
out$random_avoided <- out$random_deferred_dives/out$observed_dives
out$share_deferred <- out$deferred_dives/sum(out$deferred_dives)
out$random_share_deferred <- out$random_deferred_dives/sum(out$random_deferred_dives)
policy <- read.csv(file.path(cfg$dirs$tables,'operational_priority_cutoff.csv'))
z <- policy[policy$score=='combined',]
stopifnot(abs(1-out$avoided[out$band=='high']-z$value[z$metric=='high_retained'])<1e-12,
 abs(out$avoided[out$band=='low']-z$value[z$metric=='low_avoided'])<1e-12,
 abs(out$avoided[out$band=='ghost']-z$value[z$metric=='ghost_avoided'])<1e-12)
write_table(out,cfg,'operational_deferred_dive_breakdown.csv')

# Keep class-specific avoidance, all-dive avoidance and deferred composition
# separate: their denominators answer different management questions.
groups <- sort(unique(paste(n$reef_name,n$year,sep='||')))
g <- match(paste(d$reef_name,d$year,sep='||'),groups)
ng <- length(groups)
set.seed(cfg$hmsc$seed+27L)
weights <- cbind(rep(1,ng),replicate(500L,
 tabulate(sample.int(ng,ng,replace=TRUE),nbins=ng)))
group_sum <- function(x) as.numeric(rowsum(x,factor(g,levels=seq_len(ng)),reorder=TRUE))
ratio <- function(num,den) as.vector(crossprod(weights,num))/as.vector(crossprod(weights,den))
metrics <- list()
for (score in c('combined','random')) {
 defer <- if(score=='combined') as.numeric(d$deferred) else d$random_deferred
 all_count <- group_sum(rep(1,nrow(d)))
 deferred <- group_sum(defer)
 mid <- d$band=='intermediate'
 values <- list(all_dives_avoided=ratio(deferred,all_count),
  intermediate_retained=ratio(group_sum((1-defer)*mid),group_sum(as.numeric(mid))),
  quiet_share_deferred=ratio(group_sum(defer*(d$band %in% c('ghost','low'))),deferred))
 for(band in c('ghost','low','intermediate','high')) {
  values[[paste0('share_deferred_',band)]] <- ratio(group_sum(defer*(d$band==band)),deferred)
 }
 metrics[[score]] <- dplyr::bind_rows(lapply(names(values),function(metric) {
  v <- values[[metric]]
  data.frame(score,metric,value=v[1],low=unname(quantile(v[-1],.025)),
   high=unname(quantile(v[-1],.975)))
 }))
}
metrics <- dplyr::bind_rows(metrics)
stopifnot(abs(metrics$value[metrics$score=='combined' & metrics$metric=='all_dives_avoided']-mean(d$deferred))<1e-12)
for(score in c('combined','random')) {
 z <- metrics[metrics$score==score,]
 stopifnot(abs(sum(z$value[grepl('^share_deferred_',z$metric)])-1)<1e-12)
}
write_table(metrics,cfg,'operational_deferred_dive_metrics.csv')
print(as.data.frame(out),digits=7,row.names=FALSE)
cat('\nTotal observed dives:',nrow(d),'\nTotal deferred:',sum(d$deferred),
 '\nShare observed dives deferred:',mean(d$deferred),
 '\nRandom expected share deferred:',mean(d$random_deferred),'\n')
