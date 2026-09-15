library(Hmsc)
m <- readRDS('DataProcessing/cots_sdm_v2/iterations/hmsc_corrected_v4/models/enhanced_joint_120m.rds')
cat('XSCALEPAR_CLASS', class(m$XScalePar)[1], 'NAMES', paste(names(m$XScalePar), collapse='|'), '\n')
print(m$XScalePar)
