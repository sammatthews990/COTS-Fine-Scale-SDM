# Lizard Island clustering and Ribbon Reef isolation: network analysis plan

Status: draft, 2026-10-23. Written before the 2022/23 sub-reef connectivity
matrix (Port Douglas to Lizard) arrives, so the tests and success criteria
are fixed in advance.

## 1. Question

Wooldridge & Brodie (2015) proposed that primary COTS outbreaks start where
two conditions coincide:

- larvae survive better (nutrient-enriched water);
- the regional larval network breaks into tight local clusters in which
  nearby reefs share larvae, which happens when shelf currents are slow.

They tested this at Green Island with a 1990s reef-scale matrix (James et
al. 2002) and a binary clustering coefficient. We test the connectivity half
of the idea for the Lizard Island region with a modern, sub-reef matrix.
Stewart et al. (2026) found that dispersal network metrics alone barely
separate outbreak from non-outbreak reefs GBR-wide. Local retention and
authority score even ran opposite to expectation. A null result is therefore
plausible and worth reporting.

### Hypotheses

| | Hypothesis | Supported if | Refuted if |
|---|---|---|---|
| H1 | The Lizard group is more locally clustered than comparable reef groups elsewhere in the domain | Lizard group-retention and weighted clustering rank in the top 20 % of spatial-scan windows at a similar cross-shelf position, and the hydrodynamic matrix beats the distance-kernel control | Lizard ranks near the middle, or the distance kernel produces the same ranking |
| H2 | Clustering around Lizard is episodic: high in some releases, not others | Release-to-release variance of Lizard group-retention exceeds that of comparable windows, and high-retention releases coincide with slow local currents | Lizard's retention is stable across releases, or its peaks are unrelated to current speed |
| H3 | The Ribbon Reefs are relatively isolated from the mid-shelf | Low mid-shelf to Ribbon flux, Ribbons forming their own communities, and low in-strength at the same latitude | Ribbons are co-assigned with mid-shelf reefs and receive a mid-shelf share of supply similar to other reefs |
| H4 | Concentration hotspots in the particle-density maps match high in-strength and retention in the matrix | Spearman ρ > 0.6 between site in-strength (constant density) and day 9–10 concentration rank | Low agreement, meaning the two products need reconciling before any interpretation |

## 2. What the larval concentration maps already show

Script: `scripts/01_larval_concentration_patterns.R`. Inputs: the
`fngbr_reef_daily_means_day05/09/10` shapefiles (2,895 sub-reef sites,
113 reefs, lat −14.30 to −15.76), with 10 releases from 1 Dec 2022 to
28 Jan 2023. Comparisons use within-release percentile ranks, so differences
in total larval load between releases do not dominate.

1. **A strong north–south gradient dominates.** Reef concentration rank
   falls with distance north. Spearman ρ with latitude is −0.48 at day 5,
   −0.66 at day 9 and −0.77 at day 10 (`lc5`). Reefs north of about −14.9°,
   including the whole Lizard group, rank low at competency age.
   - This is what a domain-edge effect would produce: northern reefs lack
     upstream sources beyond the model boundary.
   - Net southward transport would produce it too.
   - **The plan must remove this before comparing Lizard with anywhere else.**
2. **The Lizard group is not a concentration hotspot.**
   - Raw ranks: it sits below other mid-shelf reefs (−0.31 at day 10).
   - Against reefs at the same latitude (latitude-detrended ranks): still
     slightly below (−0.04 to −0.05; permutation p 0.02–0.09).
   - Caveat: the detrending is partly shaped by the Lizard reefs
     themselves, so it is conservative for H1.
3. **Lizard retention is episodic.** At day 5, the Lizard group held
   1.7–1.8 times the domain mean concentration for the 28 Dec and 28 Jan
   releases. In most other releases it held well under half (`lc4`). Few
   of those peaks persist to day 9–10. This is the pattern H2 predicts, and
   it identifies which releases to examine first in the matrix.
4. **The Ribbons are relatively isolated at competency age.**
   - At the same latitude, Ribbon Reefs rank below mid-shelf reefs, and
     the gap widens with larval age: day 5 −0.04 (p 0.06), day 9 −0.05
     (p 0.01), day 10 −0.08 (p 0.0003).
   - The cross-shelf profile (`lc6`) shows why. In the central band
     (Ribbons 5–9), concentration peaks on mid-shelf reefs 20–35 km
     offshore and drops at the Ribbons.
   - In the southern band (Ribbons 1–4), Ribbons match the mid-shelf.
5. **Persistent hotspots at day 10** are the large mid-shelf reefs between
   Cooktown and Ribbon 7: Williamson, Lark, Mackay, Osterlund (`lc2`).
   Reef area does not explain this. Rank vs log area: ρ −0.10 to −0.27.

These are descriptive patterns from one season. The permutation p-values
ignore spatial autocorrelation, so they are optimistic.

## 3. Information needed from the modellers

Ask before any analysis:

- [ ] **Matrix definition.**
  - Does P_ij mean "probability that a larva released at i settles at j",
    with competency and mortality applied?
  - Which competency onset and settlement radius were used?
  - Are self-loops (P_ii) included?
- [ ] **Per-release matrices,** not only a season aggregate. Ideally the
  same 10 release dates as the concentration maps. H2 cannot be tested
  without them.
- [ ] **Source buffer.** The pilot domain ends at about −14.55°, roughly
  13 km north of Lizard Island.
  - Request releases from reefs north to at least −14.2° (Jewell, Hilder,
    Parke, Howick group), even if those reefs are not analysed.
  - Without this buffer, any Lizard result is confounded by the edge.
- [ ] **Particle mass balance per release:** the share of released particles
  that leave through the north, south and east (Coral Sea) boundaries. This
  separates edge loss from real low retention.
- [ ] **Concentration units and release design** for the concentration
  maps: are particles released uniformly per unit reef area from every
  domain reef?
- [ ] **Mean current speed and direction** at Lizard and at a central
  mid-shelf site over each release's pelagic window. These are the H2
  covariates.
- [ ] **Hydrodynamic and biological ensemble members,** if available, for
  sensitivity.
  - Choukroun et al. (2025) showed that the choice of hydrodynamic model
    changes dispersal estimates materially.

## 4. Building the networks

### Notation

- P^(r): the site-to-site settlement probability matrix for release r.
- w_i: larval production at source i. The default is constant density,
  w_i = site reef area. A sensitivity scenario uses observed COTS density
  × area, as in Stewart et al. (2026), Eq. 5.
- F_ij = w_i P_ij: settlement flux. Every metric is computed on F, so a
  large source reef is not treated like a small one.
- Season aggregate: the mean of F over releases. Report it alongside the
  per-release results, never instead of them.

### Three node scales

| Scale | Nodes | Aggregation | Role |
|---|---|---|---|
| Site | Sub-reef sites (`indexV2`; the new domain will have more) | none | Within-reef structure; links to the concentration maps and to cull/manta sites |
| Reef | GBRMPA reef (`GBRMPA_ID`) | F_RS = Σ_{i∈R} Σ_{j∈S} F_ij; P_RS = F_RS / Σ_{i∈R} w_i | Comparable with Wooldridge (2015), Kininmonth et al. (2010) and Stewart et al. (2026) |
| Reef group | Communities from the network (§5c), plus spatial-scan windows (§6) | Sum F over member sites | The direct test of "larvae retained within reef clusters" |

### Edges and thresholds

- **Weighted metrics are primary.** Binary metrics need a threshold, and
  their results depend heavily on where it is set.
- For comparison with Stewart et al. (2026) and Wooldridge, also build
  binary graphs. Use a threshold ladder: keep edges carrying the top 1, 5,
  10 and 25 % of total flux. Show each metric as a curve over the ladder
  and report only conclusions that hold across it.
- Keep self-loops for retention metrics and drop them for clustering and
  centrality.

## 5. Metrics by scale

### 5a. Site scale

| Metric | Definition | Tests |
|---|---|---|
| Local retention | P_ii | H1 (site component) |
| Own-reef supply share | Σ_{j∈reef(i)} F_ji / Σ_j F_ji | Whether supply is within-reef or between-reef |
| Constant-density in-strength | Σ_j F_ji with w = area | H4: compare with day 9–10 concentration rank |
| Weighted directed clustering | Fagiolo (2007): C_i = [(Ŵ^[1/3] + Ŵᵀ^[1/3])³]_ii / 2[d_i^tot(d_i^tot − 1) − 2d_i^↔], with Ŵ = F / max F. Report the cycle, middleman, in and out components separately | H1; the "cycle" and "in" motifs are the larval-sharing mechanism |

### 5b. Reef scale

| Metric | Definition | Tests |
|---|---|---|
| Self-retention | P_RR | H1 |
| Self-recruitment | F_RR / Σ_Q F_QR | H1; Stewart et al. found most reefs < 10 % |
| Clustering coefficient | Weighted (Fagiolo) as primary; binary over the threshold ladder for comparison with Wooldridge's CI | H1 |
| In-strength by source zone | Share of a reef's supply from inner, mid-shelf, Ribbon and own-reef sources | H3 |
| PageRank, authority | Weighted by F | Context with Stewart et al. (2026) |
| Betweenness | On −log(P), i.e. most-probable paths | Stepping-stone reefs, for example whether Lizard bridges the inner shelf and the outer barrier |
| Substantial k-step in-component | Reefs reaching R within 5 steps through edges above threshold | Comparability with Hock et al. (2014) |

### 5c. Reef-group scale

**Communities.** Use Infomap on the directed, weighted flux network (igraph
`cluster_infomap`) as the primary method. Its objective is to find where
flow stays, which matches the hypothesis.

Sensitivity checks:
- Leiden with directed modularity at several resolutions (python
  `leidenalg`; igraph's Leiden is undirected only).
- Markov stability across time scales, to show the hierarchy of clusters.

Run community detection for every release and for the season aggregate.

| Metric | Definition | Tests |
|---|---|---|
| Group retention | Σ_{i,j∈K} F_ij / Σ_{i∈K} Σ_j F_ij: the share of larvae produced in K that settle in K | H1, H2 (core metric) |
| Group self-recruitment | Σ_{i,j∈K} F_ij / Σ_i Σ_{j∈K} F_ij: the share of K's settlers produced in K | H1, H3 |
| Import share | 1 − group self-recruitment | H3: low import share means an isolated group |
| Co-assignment probability | Share of releases (and ensemble members) in which a pair of reefs falls in the same community | H2 (episodic clusters), H3 (Ribbon separation) |
| Partition stability | Adjusted Rand index / NMI between release partitions | H2 |
| Zone × zone flux matrix | F summed by inner / Lizard / mid-shelf / Ribbon / outer barrier N, row-normalised | H3, and a summary figure |
| Participation coefficient | For Ribbon nodes: how evenly their links spread across communities | H3 |

## 6. Comparison design: "compared to elsewhere"

The concentration maps show that a fixed "Lizard group" contrast is
confounded by latitude and the domain edge. The comparison therefore uses
three layers.

1. **Spatial scan (primary).**
   - Centre a window on every reef. Take all reefs within r km
     (r = 10, 15, 20, 25) and compute group retention, mean weighted
     clustering and group self-recruitment for the window.
   - Map the window values.
   - H1 is supported only if the Lizard-centred window ranks in the top
     20 % of windows at a similar distance from the mainland (±10 km), for
     at least three of the four radii.
   - This needs no hand-drawn zones.
2. **Edge buffer.**
   - Exclude windows whose centre is within 25 km of a domain boundary
     that lacks a source buffer.
   - If the modellers cannot extend sources north, report Lizard results
     as edge-affected and give them no weight in the conclusion.
   - Check that the latitude gradient seen in the concentration maps
     disappears or weakens once buffered.
3. **Null models.** Compute every metric on each null and report the
   difference from the hydrodynamic matrix.
   - **Distance kernel** (`conMatFngbrDist`, symmetric, geometry only):
     the negative control. If it gives Lizard the same clustering, the
     result reflects reef layout, not hydrodynamics.
   - **Anisotropic kernel** fitted to the hydrodynamic matrix: P as a
     function of distance and bearing. This separates "mean drift" from
     local retention structure.
   - **Flux-preserving rewiring within distance bins:** swap edges between
     pairs at similar separation, keeping each node's in- and out-strength.
     This is the null for clustering coefficients. Use 999 rewirings.

**Uncertainty.** Treat the 10 releases as replicates. Report medians and
ranges across releases, and the number of releases in which a contrast
holds. Bootstrap over releases for intervals. Do not use reef-level
permutation p-values as the main evidence, because they ignore spatial
autocorrelation.

## 7. Episodic clustering within the season (H2)

- Build a time series of Lizard window metrics across the 10 releases:
  group retention, weighted clustering, community membership. Set it
  against the distribution for comparable windows.
- Examine first the 28 Dec 2022 and 28 Jan 2023 releases, when the
  concentration maps show Lizard retention at day 5.
- **Mechanistic check.** Wooldridge & Brodie (2015) link high clustering to
  slow shelf currents. Correlate each release's Lizard group retention with
  mean current speed over the pelagic window (Spearman).
  - There are only 10 releases, so treat this as descriptive. A clear
    monotonic relationship would support H2.
- **Season context.** 2022/23 was a La Niña summer. Wooldridge & Brodie's
  Fig. 5 puts peak local clustering at neutral SOI, so this is not the
  season most likely to show the effect, and a null result in this one
  season does not refute it. Request at least one contrasting season (El
  Niño or neutral) for a fair test.

## 8. Links to other evidence

- **Concentration maps (H4).**
  - Site constant-density in-strength should track day 9–10 concentration
    rank if both products come from the same hydrodynamics. Report ρ and
    map the residuals.
  - Large disagreement means the concentration maps and the matrix use
    different biology (competency window, mortality) or releases.
- **Water quality.** The other half of Wooldridge's hypothesis is
  enriched water. Overlay summer chl-a (eReefs) on the communities and
  scan windows to find groups where high retention and high chl-a
  coincide. This is optional and descriptive.
- **COTS observations.**
  - Reef scale: relate group retention and self-recruitment to manta
    outbreak status 2–3 years after the 2022/23 spawning (2024/25–2025/26).
    With one season this is a weak, descriptive test.
  - Within-reef ranking follows the existing pilot design: the bar is a
    forward within-reef AUC of about 0.71, with the distance kernel as
    control.

## 9. Outputs and script layout

`DataProcessing/connectivity_lizard/`:

| Script | Does |
|---|---|
| `01_larval_concentration_patterns.R` | Done: concentration figures `lc1`–`lc6` and tables |
| `02_build_networks.R` | Read the matrix; pin it by SHA-256; build site, reef and zone F matrices per release and for the season |
| `03_node_metrics.R` | §5a–5b metrics, threshold ladder, null models |
| `04_communities.R` | Infomap per release, Leiden sensitivity, co-assignment, partition stability, zone × zone flux |
| `05_spatial_scan.R` | §6 scan windows, edge buffer, Lizard rank vs comparable windows |
| `06_release_variability.R` | §7 time series and current-speed check |
| `07_figures.R` | Network figures (below) |

**Planned figures:**
1. Community map per release (small multiples) plus a co-assignment
   heatmap for Lizard and Ribbon reefs.
2. Spatial-scan maps of group retention and weighted clustering, with
   Lizard's percentile marked, for the hydrodynamic matrix vs the distance
   kernel.
3. Zone × zone flux chord or matrix plot (the Ribbon isolation summary).
4. Lizard group retention by release, against the band of comparable
   windows, with current speed.
5. Clustering coefficient vs threshold curves for Lizard, other mid-shelf
   and Ribbons.
6. Concentration rank vs matrix in-strength (H4).

## 10. Caveats to state in any write-up

- **One spawning season and one hydrodynamic configuration.** Stewart &
  Bode (2025) report a median CV of about 3 for dispersal strength between
  seasons.
- **Clustering coefficients** are sensitive to thresholds, matrix
  sparsity and node scale. Compare Lizard with other parts of the same
  network, not with the absolute values reported by Wooldridge (CI about
  0.28–0.41) or Kininmonth et al. (2010).
- **Production is assumed**, not observed: constant density per reef area.
  The COTS-density scenario is a sensitivity check only.
- **Settlement is not recruitment.** Post-settlement survival and food
  (CCA, coral) are outside the network.
- **Domain edges** (§6) are the largest threat to validity for Lizard in
  particular.
