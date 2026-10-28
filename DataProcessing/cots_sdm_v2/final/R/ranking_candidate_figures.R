build_ranking_candidate_figures <- function(cfg) {
 read <- function(n) read.csv(file.path(cfg$dirs$tables,n))
 curves <- read('ranking_candidate_curves.csv')
 curves <- curves[curves$period=='2019-2026',]
 labels <- c(archived='Current ranking',mean_control='Dense-cell mean control',
  high_indicator='Prior >0.04 indicator',recent_cpue='Recency-weighted CPUE',
  hotspot_cpue='Hotspot habitat + recent CPUE',random='Random sites')
 colours <- c(archived=layer_colours[[1]],mean_control=ink$muted,
  high_indicator=layer_colours[[2]],recent_cpue=layer_colours[[3]],
  hotspot_cpue='#184f95',random=ink$primary)
 types <- c(archived='solid',mean_control='longdash',high_indicator='solid',
  recent_cpue='solid',hotspot_cpue='dotdash',random='dotted')
 plot_curve <- function(y,title,ylab,limits) {
  ggplot2::ggplot(curves,ggplot2::aes(x=.data$high_retained,y=.data[[y]],colour=.data$model,linetype=.data$model))+
   ggplot2::geom_line(linewidth=.7)+
   ggplot2::scale_colour_manual(NULL,values=colours,breaks=names(labels),labels=labels)+
   ggplot2::scale_linetype_manual(NULL,values=types,breaks=names(labels),labels=labels)+
   ggplot2::scale_x_continuous(labels=scales::label_percent(accuracy=1))+
   ggplot2::scale_y_continuous(labels=scales::label_percent(accuracy=1))+
   ggplot2::coord_cartesian(xlim=c(.75,1),ylim=limits)+
   ggplot2::labs(title=title,x='Observed >0.04 dives retained',y=ylab)+theme_report(10)
 }
 a <- plot_curve('sites_selected','A  Candidate-site coverage','Sites selected',c(.5,1))
 b <- plot_curve('dives_avoided','B  Recorded-dive reduction','Recorded dives avoided',c(0,.5))
 p <- a+b+patchwork::plot_layout(guides='collect') &
  ggplot2::theme(legend.position='bottom',legend.text=ggplot2::element_text(size=8))
 p <- p & ggplot2::guides(colour=ggplot2::guide_legend(ncol=2),linetype=ggplot2::guide_legend(ncol=2))
 save_figure(p,cfg,'fig31_ranking_candidate_tradeoffs.png',7.5,4.8)
 r <- read('ranking_candidate_assessment.csv')
 r <- r[r$rule=='habitat_only' & r$period %in% c('2019-2026','2024-2026'),]
 r$summary <- factor(r$model,levels=rev(c('h_mean','h_q50','h_q75','h_q90','h_q95')),
  labels=rev(c('Mean','Median','75th percentile','90th percentile','95th percentile')))
 r$metric <- factor(r$metric,levels=c('auc_high_all','auc_useful'),
  labels=c('>0.04 versus all other dives','Useful (>=0.02) versus low/ghost'))
 q <- ggplot2::ggplot(r,ggplot2::aes(y=.data$summary,x=.data$value,xmin=.data$low,xmax=.data$high))+
  ggplot2::geom_vline(xintercept=.5,colour=ink$muted,linetype='dotted')+
  ggplot2::geom_errorbar(orientation='y',width=.15,colour=layer_colours[[1]])+
  ggplot2::geom_point(size=2,colour=layer_colours[[1]])+
  ggplot2::facet_grid(metric~period)+ggplot2::labs(x='Within reef-year AUC',y=NULL)+theme_report(9)+
  ggplot2::theme(strip.text.y=ggplot2::element_text(angle=0),panel.spacing=grid::unit(.8,'lines'))
 save_figure(q,cfg,'fig32_habitat_polygon_percentiles.png',7.5,4.6)
}
