# 02_network_starter.R
#
# Question: Is the Lizard Island area more tightly "clustered" (larvae shared
#   and kept among nearby reefs) than elsewhere, and are the Ribbon Reefs cut
#   off from the mid-shelf?
# Method: Starter version for the RA (see RA_BRIEF.md). Builds a reef-scale
#   larval flux network from a sub-reef connectivity matrix and computes:
#   self-retention, self-recruitment, weighted clustering, Infomap
#   communities, a zone-to-zone flux table and a 20 km moving-window
#   "group retention" scan.
# Validation: None yet. Currently runs on conMatFngbrDist (a distance-only
#   kernel, NOT hydrodynamics), so results are a dry run and the future
#   negative control. Swap matrix_file when the 2022/23 hydro matrix arrives.
#
# Run from the repository root (after script 01, which writes reef_zones.csv):
#   Rscript DataProcessing/connectivity_lizard/scripts/02_network_starter.R

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(igraph)
})

root <- "DataProcessing/connectivity_lizard"
dir_tab <- file.path(root, "outputs", "tables")
matrix_file <- "DataProcessing/data/conMatFngbrDist 1.xlsx"   # <- swap for the hydro matrix
matrix_label <- "distance_kernel"                              # <- e.g. "hydro_2022_23"
sites_file <- "DataProcessing/data/fngbr_LarvalConcentrations/fngbr_reef_daily_means_day10.shp"
edge_keep <- 0.10      # keep edges carrying the top 10 % of flux for graph metrics
scan_radius_km <- 20
lizard_reef <- "14-116d" # Lizard Island Reef (Lagoon)

out <- function(x, name) {
  path <- file.path(dir_tab, paste0("net_", matrix_label, "_", name, ".csv"))
  write.csv(x, path, row.names = FALSE)
  message("Wrote ", path)
}

# ---- 1. Matrix ---------------------------------------------------------------
# Convention used here: rows = source site, columns = destination site,
# P[i, j] = probability a larva released at i settles at j.
# CHECK THIS with the modellers for every new matrix.
P <- as.matrix(readxl::read_excel(matrix_file, sheet = 1, col_names = FALSE,
                                  .name_repair = "minimal"))
dimnames(P) <- NULL
stopifnot(nrow(P) == ncol(P))
message("Matrix ", nrow(P), " x ", ncol(P), "; symmetric: ", isSymmetric(P),
        "; share of zero entries: ", round(mean(P == 0), 3))

# ---- 2. Sites -> reefs ---------------------------------------------------------
# Matrix row/col k is sub-reef site indexV2 == k.
sites <- st_read(sites_file, quiet = TRUE) |>
  st_drop_geometry() |>
  group_by(indexV2) |>
  summarise(GBRMPA_ID = first(GBRMPA_ID), area = sum(area), .groups = "drop") |>
  arrange(indexV2)
stopifnot(identical(sites$indexV2, seq_len(nrow(P))))

zones <- read.csv(file.path(dir_tab, "reef_zones.csv"))  # from script 01
reefs <- zones |> select(GBRMPA_ID, name, zone, lon, lat) |>
  filter(GBRMPA_ID %in% sites$GBRMPA_ID)

# Flux: production at the source (constant density => proportional to area)
# times settlement probability.
Fsite <- sites$area * P

Z <- outer(sites$GBRMPA_ID, reefs$GBRMPA_ID, `==`) * 1   # site x reef membership
Freef <- t(Z) %*% Fsite %*% Z                              # reef x reef flux
dimnames(Freef) <- list(reefs$GBRMPA_ID, reefs$GBRMPA_ID)
production <- as.vector(t(Z) %*% sites$area)

# ---- 3. Reef metrics -----------------------------------------------------------
# Weighted directed clustering (Fagiolo 2007, "total" version): how strongly
# a reef's neighbours are themselves linked, using edge weights.
fagiolo_clustering <- function(W) {
  diag(W) <- 0
  A <- (W > 0) * 1
  W3 <- (W / max(W))^(1 / 3)
  S <- W3 + t(W3)
  num <- diag(S %*% S %*% S)
  d_tot <- rowSums(A) + colSums(A)
  d_bi <- diag(A %*% A)
  num / (2 * (d_tot * (d_tot - 1) - 2 * d_bi))
}

# Keep only the strongest edges (top edge_keep share of between-reef flux).
Fb <- Freef; diag(Fb) <- 0
cut <- quantile(Fb[Fb > 0], 1 - edge_keep)
Fthin <- Fb * (Fb >= cut)

reef_metrics <- reefs |>
  mutate(
    self_retention = diag(Freef) / production,           # share of own larvae kept
    self_recruitment = diag(Freef) / colSums(Freef),    # share of settlers home-grown
    in_strength = colSums(Fb) / production,              # imported supply per unit area
    clustering = fagiolo_clustering(Fthin)
  )

# ---- 4. Communities ------------------------------------------------------------
g <- graph_from_adjacency_matrix(Fthin, mode = "directed", weighted = TRUE)
set.seed(1)
comm <- cluster_infomap(g, e.weights = E(g)$weight)
reef_metrics$community <- membership(comm)[reef_metrics$GBRMPA_ID]
message("Infomap found ", length(unique(reef_metrics$community)), " communities")
out(reef_metrics, "reef_metrics")

# ---- 5. Zone-to-zone flux (Ribbon isolation) -----------------------------------
Zz <- outer(reefs$zone, sort(unique(reefs$zone)), `==`) * 1
Fzone <- t(Zz) %*% Fb %*% Zz
dimnames(Fzone) <- list(from = sort(unique(reefs$zone)), to = sort(unique(reefs$zone)))
# Column-normalised: where does each zone's incoming supply come from?
zone_supply <- as.data.frame(as.table(sweep(Fzone, 2, colSums(Fzone), "/")),
                             responseName = "share_of_to_zone_supply")
out(zone_supply, "zone_supply")

# ---- 6. Moving-window scan ---------------------------------------------------------
# For each reef, take every reef within scan_radius_km and ask: of the larvae
# produced in that window, what share settles back inside it?
xy <- st_coordinates(st_transform(st_as_sf(reefs, coords = c("lon", "lat"), crs = 4326), 32755))
D <- as.matrix(dist(xy)) / 1000
group_retention <- function(m) sum(Freef[m, m]) / sum(Freef[m, ])
scan <- reefs |>
  mutate(n_reefs = colSums(D <= scan_radius_km),
         group_retention = vapply(seq_len(nrow(reefs)),
                                  \(k) group_retention(which(D[k, ] <= scan_radius_km)), 0),
         pct_rank = percent_rank(group_retention))
out(scan, "scan_20km")
liz <- scan[scan$GBRMPA_ID == lizard_reef, ]
message(sprintf("Lizard-centred %d km window: group retention %.3f, percentile %.0f of %d windows",
                scan_radius_km, liz$group_retention, 100 * liz$pct_rank, nrow(scan)))
