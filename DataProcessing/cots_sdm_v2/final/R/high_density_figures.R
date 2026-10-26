# Trade-off for the unchanged forward ranking, evaluated at >0.04 CPUE.
build_high_density_figures <- function(cfg) {
 d <- utils::read.csv(file.path(cfg$dirs$tables,'cull_high_density_fixed.csv'))
 d <- d[d$design=='forward' & d$period=='2019-2026' & d$score %in% c('habitat','combined','random'),]
 x <- d[d$metric=='sites_searched',c('score','nominal_fraction','value')]
 names(x)[3] <- 'sites'
 d <- merge(d[d$metric=='high_retained',],x,by=c('score','nominal_fraction'))
 labels <- c(habitat='Habitat alone',combined='Habitat + recent evidence',random='Random, same site effort')
 d$ranking <- factor(labels[d$score],levels=labels)
 colours <- stats::setNames(c(layer_colours[[1]],layer_colours[[3]],ink$secondary),labels)
 p <- ggplot2::ggplot(d,ggplot2::aes(sites,value,colour=ranking,linetype=ranking))+
  ggplot2::geom_ribbon(ggplot2::aes(ymin=low,ymax=high,fill=ranking),alpha=.1,colour=NA,linetype=0)+
  ggplot2::geom_line(linewidth=.8)+
  ggplot2::geom_hline(yintercept=c(.8,.9,.95,.99),colour=ink$muted,linetype='dotted',linewidth=.35)+
  ggplot2::scale_colour_manual(values=colours)+ggplot2::scale_fill_manual(values=colours)+
  ggplot2::scale_linetype_manual(values=c('solid','longdash','dotted'))+
  ggplot2::scale_x_continuous(labels=scales::label_percent(),limits=c(0,1))+
  ggplot2::scale_y_continuous(labels=scales::label_percent(),limits=c(0,1))+
  ggplot2::labs(x='Designated site-years searched',y='Dives above 0.04 CPUE retained',colour=NULL,fill=NULL,linetype=NULL,
    title='Search reductions depend on high-density retention',
    subtitle='Existing forward evidence ranking, 2019–2026; no model refit')+
  theme_report(10)+ggplot2::theme(legend.position='bottom')
 ggplot2::ggsave(file.path(cfg$dirs$figures,'fig20_high_density_retention.png'),p,width=7.5,height=5,dpi=300,bg=ink$surface)
}
