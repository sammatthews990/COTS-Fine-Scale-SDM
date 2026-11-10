# Reader-facing surveillance summary from the existing frozen aggregation table.
# Original detailed comparisons and figure 13 remain in the appendices.
build_report_highlight_figures <- function(cfg) {
 d <- utils::read.csv(file.path(cfg$dirs$tables,'manta_aggregation_cv_summary.csv'))
 d <- d[d$subset=='aggregations (>= 2 strong tows), all years',]
 labels <- c(layer_cv='Habitat priority',systematic='Evenly spaced tows',random='Random selection')
 long <- dplyr::bind_rows(lapply(names(labels),function(n) data.frame(
  retained=d$keep, found=d[[n]], low=d[[paste0(n,'_low')]], high=d[[paste0(n,'_high')]],
  strategy=labels[[n]])))
 colours <- stats::setNames(c(layer_colours[[2]],ink$secondary,ink$muted),labels)
 long$strategy <- factor(long$strategy,levels=labels)
 p <- ggplot2::ggplot(long,ggplot2::aes(retained,found,colour=strategy,linetype=strategy))+
  ggplot2::geom_ribbon(ggplot2::aes(ymin=low,ymax=high,fill=strategy),alpha=.1,colour=NA,linetype=0)+
  ggplot2::geom_line(linewidth=.85)+
  ggplot2::geom_hline(yintercept=.95,colour=ink$muted,linetype='dotted',linewidth=.35)+
  ggplot2::scale_colour_manual(values=colours)+ggplot2::scale_fill_manual(values=colours)+
  ggplot2::scale_linetype_manual(values=c('longdash','solid','dotted'))+
  ggplot2::scale_x_continuous(labels=scales::label_percent(),limits=c(.3,1))+
  ggplot2::scale_y_continuous(labels=scales::label_percent(),limits=c(.4,1))+
  ggplot2::labs(x='Tows retained in each survey',y='Observed aggregations still detected',
    colour=NULL,fill=NULL,linetype=NULL,title='Evenly spaced tows retain more aggregations')+
  theme_report(10)+ggplot2::theme(legend.position='bottom')
 ggplot2::ggsave(file.path(cfg$dirs$figures,'fig21_manta_surveillance_highlight.png'),p,width=7.5,height=4.8,dpi=300,bg=ink$surface)
}
