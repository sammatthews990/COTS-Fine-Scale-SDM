# Geomorphology-led Hmsc iteration 03

This directory contains an independent successor experiment. It reads shared
source data and frozen v2 artifacts, but it does not modify the legacy model or
use the legacy linear predictor, probabilities, or ranks during fitting.

The iteration compares habitat models at their fitted support (60 m, 120 m,
and track-footprint weighted) and keeps reef-year outbreak probability as a
separate activation product. Legacy clean + year is used only in benchmark
tables and expert comparison panels.

Run entry points with `Rscript scripts/<script>.R` from the repository root.

