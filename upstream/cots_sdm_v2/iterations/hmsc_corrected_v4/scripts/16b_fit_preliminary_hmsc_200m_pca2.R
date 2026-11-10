args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
v4_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
source_file <- file.path(v4_dir, "scripts", "04_fit_preliminary_hmsc_120m.R")
code <- readLines(source_file, warn = FALSE)
code <- gsub("120m", "200m", code, fixed = TRUE)
start <- grep("^variants <- list\\($", code)
end <- grep("^\\)$", code)
end <- end[end > start][1]
replacement <- c(
  'variants <- list(',
  '  enhanced_cots_200support_pca2 = list(',
  '    domain = "enhanced", responses = cots_responses,',
  '    base_terms = c("geomorphic_class", "AhyaD", "Aspat", "Aten", "turbid", "waves_Hs", "waves_Tp")',
  '  )',
  ')'
)
code <- c(code[seq_len(start - 1L)], replacement, code[(end + 1L):length(code)])
code <- gsub(
  "fit_family_pca(data, variance = 0, max_components = 1L)",
  "fit_family_pca(data, variance = 0, max_components = 2L)",
  code, fixed = TRUE
)
code <- gsub("support_m = 120L", "support_m = 200L", code, fixed = TRUE)
eval(parse(text = code, keep.source = TRUE), envir = globalenv())
