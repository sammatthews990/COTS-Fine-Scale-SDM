library(Hmsc)
m <- readRDS('DataProcessing/cots_sdm_v2/iterations/hmsc_corrected_v4/models/enhanced_joint_120m.rds')
cat('NAMES', paste(names(m), collapse='|'), '\n')
for (nm in intersect(c('X','XScaled','Tr','Y','YScaled','XData','XFormula'), names(m))) {
  z <- m[[nm]]
  cat(nm, class(z)[1], paste(dim(z), collapse='x'), '\n')
}
cat('scaleX', m$scaleX, 'XScale', paste(names(m$XScale), collapse='|'), '\n')
if (!is.null(m$XScale)) print(m$XScale)
