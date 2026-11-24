# Observation footprint and northern prototype for the supplied report template.
# This is a location map, not a new validation or model fit.
build_ccip_report_figures <- function(cfg) {
 d <- readRDS(checked_input(cfg,'assessment_dives'))
 d <- d[is.finite(d$longitude) & is.finite(d$latitude) & is.finite(d$total_cots) &
  is.finite(d$bottom_time) & d$bottom_time>0 & d$longitude>140 & d$longitude<155 & d$latitude< -8 & d$latitude> -27,]
 # One symbol per observed reef avoids implying a probability surface.
 pts <- d |> dplyr::group_by(.data$reef_name) |>
  dplyr::summarise(lon=stats::median(.data$longitude),lat=stats::median(.data$latitude),.groups='drop')
 path <- checked_input(cfg,'reef_features')
 layer <- sf::st_layers(path)$name[1]
 land <- sf::st_read(path,query=paste0('SELECT * FROM "',layer,'" WHERE FEAT_NAME IN (\'Mainland\',\'Island\')'),quiet=TRUE)
 land <- sf::st_transform(land,4326)
 land <- suppressWarnings(sf::st_crop(land,c(xmin=140,xmax=155,ymin=-27,ymax=-8)))
 reefs <- read_designated_sites(cfg)
 cr <- reefs[grepl('Lizard Island|Eyrie|MacGillivray|North Direction',reefs$reef_name),]
 cr <- sf::st_transform(cr,4326)
 lab <- pts[grepl('Lizard Island|Eyrie|MacGillivray|North Direction',pts$reef_name),]
 lab$reef <- ifelse(grepl('Lizard',lab$reef_name),'Lizard Island',ifelse(grepl('Eyrie',lab$reef_name),'Eyrie',ifelse(grepl('MacGillivray',lab$reef_name),'MacGillivray','North Direction')))
 lab <- lab |> dplyr::group_by(.data$reef) |> dplyr::summarise(lon=mean(.data$lon),lat=mean(.data$lat),.groups='drop')
 left <- ggplot2::ggplot()+ggplot2::geom_sf(data=land,fill=ink$land,colour=NA)+
  ggplot2::geom_point(data=pts,ggplot2::aes(lon,lat),colour=layer_colours[[1]],size=.55,alpha=.7)+
  ggplot2::annotate('rect',xmin=145.2,xmax=145.65,ymin=-14.9,ymax=-14.5,fill=NA,colour=ink$primary,linewidth=.6)+
  ggplot2::coord_sf(xlim=c(141,153.5),ylim=c(-25,-9.5),expand=FALSE)+
  ggplot2::labs(tag='a',title='Cull observation footprint',x=NULL,y=NULL)+theme_report(10)
 right <- ggplot2::ggplot()+ggplot2::geom_sf(data=land,fill=ink$land,colour=NA)+
  ggplot2::geom_sf(data=cr,fill=layer_colours[[1]],colour=NA,alpha=.4)+
  ggplot2::geom_point(data=lab,ggplot2::aes(lon,lat),size=1.5,colour=ink$primary)+
  ggplot2::geom_text(data=lab,ggplot2::aes(lon,lat,label=reef),nudge_x=.018,hjust=0,size=3.2)+
  ggplot2::scale_x_continuous(breaks=c(145.3,145.45,145.6))+
  ggplot2::scale_y_continuous(breaks=c(-14.5,-14.6,-14.7,-14.8))+
  ggplot2::coord_sf(xlim=c(145.25,145.65),ylim=c(-14.85,-14.5),expand=FALSE)+
  ggplot2::labs(tag='b',title='Northern calibration reefs',x=NULL,y=NULL)+theme_report(10)
 domain <- sf::st_read(checked_input(cfg,'report_metapop_domain'),quiet=TRUE)
 stopifnot(nrow(domain)==2895,length(unique(domain$GBRMPA_ID))==113)
 domain <- sf::st_transform(domain,4326)
 # Panel c: the COMPAS-EMS model mesh (whole domain, with the mesh around
 # Lizard Island as an inset), supplied by the author on 20 November 2026 as
 # images in report/assets.
 mesh_dir <- file.path(cfg$final_dir, 'report', 'assets')
 img_panel <- function(f) {
  x <- png::readPNG(file.path(mesh_dir, f))
  ggplot2::ggplot()+ggplot2::annotation_raster(x, 0, ncol(x), 0, nrow(x))+
   ggplot2::coord_fixed(xlim=c(0, ncol(x)), ylim=c(0, nrow(x)), expand=FALSE)+theme_map(10)
 }
 mesh_zoom <- img_panel('hydro_mesh_lizard.png')+
  ggplot2::theme(plot.background=ggplot2::element_rect(fill=ink$surface, colour=ink$primary, linewidth=.5),
   plot.margin=ggplot2::margin(1, 1, 1, 1))
 pending <- img_panel('hydro_mesh_domain.png')+ggplot2::labs(tag='c', title='Hydrodynamic model mesh')+
  patchwork::inset_element(mesh_zoom, left=0, bottom=0, right=.56, top=.28, align_to='panel')
 full_domain <- ggplot2::ggplot()+
  ggplot2::geom_sf(data=land,fill=ink$land,colour=NA)+
  ggplot2::geom_sf(data=domain,fill=layer_colours[[1]],colour=ink$primary,linewidth=.12)+
  ggplot2::geom_point(data=lab,ggplot2::aes(lon,lat),size=1,colour=ink$primary)+
  ggplot2::coord_sf(xlim=c(145.05,145.9),ylim=c(-15.82,-14.25),expand=FALSE)+
  ggplot2::labs(tag='d',title='Population-model domain',x=NULL,y=NULL)+theme_report(10)
 p <- (left+right)/(pending+full_domain)+patchwork::plot_layout(widths=c(1,1.1))
 ggplot2::ggsave(file.path(cfg$dirs$figures,'fig29_report_study_area.png'),p,width=7.5,height=8,dpi=300,bg=ink$surface)
}
