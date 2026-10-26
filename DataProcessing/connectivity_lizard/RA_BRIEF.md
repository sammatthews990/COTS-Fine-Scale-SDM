# Lizard Island larval network: RA brief

This is the informal, cut-down version of [ANALYSIS_PLAN.md](ANALYSIS_PLAN.md).
Read this first; dip into the plan when you want the detail behind a step.

## The idea in a nutshell

There is an old hypothesis (Wooldridge & Brodie 2015) about how
crown-of-thorns starfish outbreaks get started. Most years, larvae get
flushed along the reef. In some years currents are slow, and groups of
nearby reefs end up trading larvae back and forth among themselves. If that
"huddle" also sits in food-rich water, enough larvae survive to kick off an
outbreak. They tested this at Green Island with 1990s models.

We want to test the connectivity half for the **Lizard Island** area,
using a new, much finer connectivity model for the 2022/23 spawning season.
There are two questions:

1. **Is the Lizard area more of a "huddle"** than other parts of the region?
   That is, do reefs near Lizard keep and share more of their larvae than
   reef groups elsewhere?
2. **Are the Ribbon Reefs (the outer barrier) cut off** from the mid-shelf
   reefs, getting few larvae from them?

A "no" to either is a perfectly good result. Another recent paper (Stewart
et al. 2026) found network metrics on their own explain
outbreaks poorly. We just need to answer the questions cleanly.

## What we already know

Script `01_larval_concentration_patterns.R` used maps of modelled larval
concentration (how many larvae are sitting over each patch of reef, 5, 9
and 10 days after release, for 10 release dates in Dec 2022–Jan 2023). The
figures are in `outputs/figures/` (`lc1`–`lc6`). Headlines:

- **North is low, south is high.** Larval concentration drops steadily as
  you go north (`lc5`). Some of this is probably an artefact: the model
  domain ends just north of Lizard, so reefs up there have no "upstream"
  reefs feeding them. Keep this in mind for everything below.
- **Lizard doesn't stand out**, even when compared only with reefs at the
  same latitude.
- **But Lizard has occasional moments.** For the 28 Dec and 28 Jan
  releases, the Lizard reefs held far more larvae than average at day 5
  (`lc4`). This is the kind of "sometimes it huddles" pattern the
  hypothesis predicts.
- **The Ribbons get fewer larvae than mid-shelf reefs at the same
  latitude,** especially by day 10 (`lc6`, middle row). So there's an early
  hint that they are somewhat isolated.

## Getting set up

Everything runs in R from the repository root:

```
Rscript DataProcessing/connectivity_lizard/scripts/01_larval_concentration_patterns.R
Rscript DataProcessing/connectivity_lizard/scripts/02_network_starter.R
```

Script 01 makes the figures and `outputs/tables/reef_zones.csv`. This file
puts each of the 113 reefs in a zone: inner shelf, Lizard group, other
mid-shelf, Ribbon Reefs, or outer barrier north of Lizard.

Script 02 is your starting point for the network work. Right now it runs on
`conMatFngbrDist 1.xlsx`. That is **not** a real connectivity matrix: it
just says "nearby reefs are connected", based on distance alone. That is
deliberate.
- It lets you get the code working before the real matrix arrives.
- It later becomes our **control**. If the real model only shows what
  distance alone shows, the hydrodynamics aren't adding anything.

When the real matrix arrives you change two lines at the top of 02
(`matrix_file` and `matrix_label`) and run it again.

Packages: `sf`, `dplyr`, `igraph`, `readxl`, `ggplot2` (all installed).

## What script 02 does

1. **Reads the matrix.**
   - Rows are source sites and columns are destination sites. Each of the
     2,895 rows/columns is a sub-reef site (`indexV2`).
   - `P[i, j]` = the chance a larva from site i settles at site j.
2. **Turns probabilities into larval "flux".**
   - It multiplies each row by the source site's area, on the assumption
     that every reef patch releases the same number of larvae per m².
   - It then adds sites up into reefs, giving a 113 × 113 reef-to-reef
     flux matrix.
3. **Works out per-reef numbers** (`net_<label>_reef_metrics.csv`):

   | Column | Plain English |
   |---|---|
   | `self_retention` | Of the larvae this reef releases, what share come back to it? |
   | `self_recruitment` | Of the larvae settling here, what share were home-grown? |
   | `in_strength` | How many larvae arrive from *other* reefs, per unit of this reef's area? |
   | `clustering` | Do this reef's strong partners also trade strongly with each other? High = part of a tight huddle. `NA` = no strong links at all. |
   | `community` | Which "neighbourhood" of reefs it belongs to. Infomap finds groups where larvae tend to stay. |

4. **Zone-to-zone table** (`net_<label>_zone_supply.csv`): for each zone,
   where its incoming larvae come from. This is the main Ribbon-isolation
   number.
5. **Moving-window scan** (`net_<label>_scan_20km.csv`):
   - Around every reef, it takes all reefs within 20 km and asks: of the
     larvae released in this circle, what share settle back inside it?
     That share is the "group retention".
   - It then ranks every circle. The Lizard-centred circle's rank is the
     main answer to question 1.
   - Dry run on the distance-only matrix: Lizard comes out at the 13th
     percentile. So reef layout alone does *not* make Lizard look
     huddled. The real model has to show something distance can't.

## Your tasks, roughly in order

### Task 1: get comfortable (now, with the distance-only matrix)

- Run 01 and 02. Look at the figures and open the CSVs.
- **Make a map of the scan results:** colour each reef by `group_retention`
  and mark Lizard. Copy the map style from script 01 (`base_map()`,
  `ramp_blue`).
- **Make a map of the Infomap communities.**
- **Try other window sizes** (10, 15 and 25 km; it's one variable,
  `scan_radius_km`). Does Lizard's rank move around a lot?
- **Try other edge thresholds** (`edge_keep` = 0.05, 0.25). The clustering
  numbers depend heavily on this, which is why we lean on group retention
  instead.

### Task 2: sanity-check the real matrix when it arrives

Before computing anything, check these. Write down what you find.

- **Orientation.** Is it really rows = source? Ask the modellers; don't
  guess. For a real hydro model, `P[i, j]` and `P[j, i]` should differ.
  `isSymmetric(P)` should be `FALSE`. If it's `TRUE`, something is off.
- **Row sums** should be ≤ 1, because a larva can't settle more than once.
  Many will be well below 1: larvae die or drift away.
- **Indexing.** Same site numbering as our shapefile? The new domain
  (Port Douglas to Lizard) is bigger, so you'll need its site shapefile.
- **Per-release or season total?** We want one matrix per release date if
  possible.
- **Does it agree with the concentration maps?**
  - Compute each site's incoming flux, `colSums(Fsite)`.
  - Correlate it (Spearman) with the day 9–10 concentration rank from
    script 01.
  - If they disagree badly, flag it before going further.

### Task 3: question 1, is Lizard a huddle?

- Re-run 02 on the real matrix. Rerun the scan at 10, 15, 20 and 25 km.
- Compare Lizard's window **with windows at a similar distance from the
  coast**, not with everything. Inner-shelf and outer-barrier circles
  behave differently. `dist_coast_km` is in `reef_zones.csv`.
- **Drop circles that touch the edge of the model domain** (within about
  25 km). They're missing incoming larvae. If Lizard itself is too close
  to the edge, say so. That is a real finding about the data, not a
  failure.
- Put the real-matrix result next to the distance-only result. The
  interesting part is the *difference*.
- A useful one-line summary: "Lizard ranks X out of Y comparable windows
  for group retention (distance-only: Z)".

### Task 4: question 2, are the Ribbons cut off?

- **Zone supply table.** What share of the Ribbons' incoming larvae comes
  from mid-shelf reefs, compared with other zones? Compare with the
  distance-only version.
- **Communities.** Do Ribbon reefs form their own communities, or get
  lumped with mid-shelf reefs?
- **Split by latitude.** The concentration maps say the answer differs
  north and south: Ribbons 1–4 behave differently from Ribbons 5–9. The
  latitude bands are in script 01 (`band_breaks`).

### Task 5 (if we get per-release matrices): does Lizard huddle only sometimes?

- Loop the scan over each release's matrix and plot Lizard's group
  retention by release date. Show the range of comparable windows as a
  grey band behind it.
- Look first at **28 Dec and 28 Jan**, the releases where the
  concentration maps showed larvae piling up at Lizard.
- If the modellers can give us average current speed per release, plot
  Lizard's retention against it. The hypothesis says slow currents mean
  more huddling. With only 10 points this is a "look at it" plot, not a
  statistics test.

## Things that will trip you up

- **Edge effects.** The biggest one. Reefs near the edge of the model look
  isolated because their neighbours over the boundary don't exist in the
  model.
- **Distance does a lot of the work.** Nearby reefs are always more
  connected. Always ask "is this more than distance alone would give?"
- **Thresholds.** Anything that needs a cut-off for "strong link"
  (clustering, communities) can change a lot with the cut-off. Try a few
  values and only trust what holds across them.
- **Raw probabilities vs flux.** Big reefs release more larvae. Use the
  area-weighted flux (`Freef`), not `P`, unless you have a reason not to.
- **One season only, and a La Niña one.** The hypothesis says huddling is
  strongest in neutral years, so a "no" from 2022/23 doesn't kill it.
- **p-values.** Reefs next to each other aren't independent, so classic
  p-values will look better than they should. Lean on ranks against
  comparable windows, the distance-only comparison, and consistency across
  release dates.

## What to send back

- [ ] A short notes file (half a page is fine): what you ran, the sanity
  checks, the headline numbers, anything weird.
- [ ] Maps: scan group retention (real vs distance-only) and communities.
- [ ] The zone supply table, real vs distance-only.
- [ ] A Lizard-by-release plot, if we have per-release matrices.
- [ ] Code committed to a branch, with any new scripts numbered `03_...`
  and so on. Give each script a short header saying what question it
  answers, like 01 and 02.

## Mini-glossary

- **Connectivity matrix:** a big table of "chance a larva from here ends up
  there", for every pair of sites.
- **Flux:** larvae released × that chance. It's the "how many larvae move
  from A to B" version.
- **Self-retention / local retention:** larvae that come home to the reef
  they were released from.
- **Group retention:** the same idea for a group of reefs. Of the larvae
  the group releases, what share stays inside the group?
- **Clustering coefficient:** whether your partners are also each other's
  partners. High = tight-knit group.
- **Community (Infomap):** a group of reefs where larvae tend to circulate
  among themselves more than they leave.
- **Distance kernel:** a fake connectivity matrix based only on how far
  apart reefs are. It's our baseline.

## Reading (ask Sam for the PDFs)

- Wooldridge & Brodie 2015, *Mar Pollut Bull*. The hypothesis; Figs 4–5
  are the key ones.
- Stewart et al. 2026, *Coral Reefs*. The modern GBR-wide test. Useful for
  metric definitions (their Eqs 1–5).
