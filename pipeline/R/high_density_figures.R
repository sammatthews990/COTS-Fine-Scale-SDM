# Trade-off for the unchanged forward ranking, evaluated at >0.04 CPUE.
build_high_density_figures <- function(cfg) {
 d <- utils::read.csv(file.path(cfg$dirs$tables,'cull_high_density_fixed.csv'))
 d <- d[d$design=='forward' & d$period=='2019-2026' & d$score %in% c('habitat','combined','random'),]
 x <- d[d$metric=='sites_searched',c('score','nominal_fraction','value')]
 names(x)[3] <- 'sites'
 d <- merge(d[d$metric=='high_retained',],x,by=c('score','nominal_fraction'))
 labels <- c(habitat='Habitat alone',combined='Habitat + recent evidence',random='Random')
 d$ranking <- factor(labels[d$score],levels=labels)
 colours <- stats::setNames(unname(key_colours[c('habitat','combined','random')]),labels)
 p <- ggplot2::ggplot(d,ggplot2::aes(sites,value,colour=ranking,linetype=ranking))+
  ggplot2::geom_ribbon(ggplot2::aes(ymin=low,ymax=high,fill=ranking),alpha=.1,colour=NA,linetype=0)+
  ggplot2::geom_line(linewidth=.8)+
  ggplot2::geom_hline(yintercept=c(.8,.9,.95,.99),colour=ink$muted,linetype='dotted',linewidth=.35)+
  ggplot2::scale_colour_manual(values=colours)+ggplot2::scale_fill_manual(values=colours)+
  ggplot2::scale_linetype_manual(values=c('solid','solid','dashed'))+
  ggplot2::scale_x_continuous(labels=scales::label_percent(),limits=c(0,1))+
  ggplot2::scale_y_continuous(labels=scales::label_percent(),limits=c(0,1))+
  ggplot2::labs(x='Designated site-years searched',y='High-density dives kept',colour=NULL,fill=NULL,linetype=NULL)+
  theme_report(9)+ggplot2::theme(legend.position='inside',legend.position.inside=c(.98,.03),legend.justification=c(1,0),
    legend.key.width=grid::unit(.8,'cm'))
 # The product's integer A-C rule differs slightly from the percentile grid
 # used for matched-retention interpolation. Mark its exact measured result.
 op <- utils::read.csv(file.path(cfg$dirs$tables,'operational_priority_cutoff.csv'))
 op <- op[op$score=='combined',]
 xop <- op$value[op$metric=='sites_selected']; yop <- op$value[op$metric=='high_retained']
 p <- p+ggplot2::annotate('point',x=xop,y=yop,shape=23,size=3,fill=ink$surface,colour=ink$primary)
 save_figure(p,cfg,'fig20_high_density_retention.png',width=5.5,height=4.2)
}
