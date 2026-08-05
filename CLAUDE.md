# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

An R analysis codebase (not a package) for modelling the spread and occupancy of Red Dwarf Honeybees (RDHB), an invasive species, using surveillance survey data from an eradication response. There is no build system, package manifest, or test suite — this is a collection of sequential/independent R scripts run interactively or via `Rscript`, plus R Markdown reports.

## Data access

Survey data are not in the repo (not publicly available). Scripts expect a `DATA_PATH` env var pointing at the local data directory, set via:

```r
usethis::edit_r_environ()
# then add: DATA_PATH="your/local/path/to/data/directory"
```

Restart R after editing. `src/a-setup.R` builds `data_dir <- file.path(Sys.getenv("DATA_PATH"), "RDHB")`, which downstream scripts read from. Without access to this data directory, scripts that call `source("src/b-data-organisation.R")` (directly or transitively) cannot be run end-to-end.

## Running scripts

There is no CLI entry point — open `RDHB_Occupancy.Rproj` in RStudio, or run individual scripts from the repo root with `Rscript src/<script>.R` (working directory must be the repo root; scripts use relative paths like `"src/a-setup.R"` and `"out/..."`).

`src/z-monthly-update.R` is the closest thing to a "main" pipeline: it sets `current.data` (the surveillance `.xlsx` filename to use for that run), sources data organisation, regenerates the CpUE figure and the multi-season occupancy model/figures, and re-renders `ms/RDHB_report.Rmd` to a dated HTML file. Editing `current.data` at the top of this script is how you point analyses at a new data export.

## Architecture: three independent modelling lineages

All scripts live flat in `/src`, prefixed with letters that indicate both order and lineage. Every modelling script starts by sourcing `src/a-setup.R` (shared helpers/config) and `src/b-data-organisation.R` (loads and cleans the raw surveillance `.xlsx` into a point-level `sf` dataframe `df`, plus a separate colony-level `cny.df`). From there, the three lineages diverge and do not depend on each other:

- **`a`–`f`: occupancy models of increasing complexity**, fit with `unmarked`/`spOccupancy` on data aggregated into a spatial grid and time steps via `aggregate_data()` (from `a-setup.R`). `c` = single-season, `d` = time-varying detection, `e` = dynamic (multi-season) occupancy, `f` = multi-season occupancy with spatial + temporal autocorrelation (via `spOccupancy::stPGOcc`). **`f` is the most actively used/useful model of this group** and is what `z-monthly-update.R` runs. `f-1-...data-organisation.R` is `f`'s dedicated data-prep step.

- **`g`–`h`: dynamic spatial occupancy model**, an experimental extension meant to fold control effort (colony/hive removal) into the occupancy framework via distance-weighted colonisation/extinction dynamics. `g-1` defines the shared process-model functions (colonisation/extinction as a function of distance-weighted neighbour occupancy), `g-2`/`g-3`/`g-4` simulate data and test-fit the model in `greta` and JAGS respectively before touching real data, `h-1` organises real data into the required format, `h-...JAGS.R` / `h-...Nimble.R` fit the model itself via two alternate backends. This lineage is described in the README as "kind of works, but needs more development" — treat it as exploratory, not production.

- **`i`–`j`: point-process ("colony-location") model**, a different modelling paradigm that infers the *locations* of colonies from imperfect detections of foraging bees plus a foraging-distance kernel, rather than modelling grid-cell occupancy directly. `i-0` has shared simulation functions, `i-1` simulates fake data, `i-2-pp_*` fits the simulated case in Nimble/JAGS to validate the approach. `j-pp-static-get-data(.R/_pawsey.R)` prepares real survey + foraging data (transformed to UTM 50S coordinates, filtered to a rolling time window via `up.to.date`), `j-pp-static.R` fits the real static (single time-window) model in JAGS against `src/model-files/pp-static-JAGS*.txt`, and `colonies_compare.R` validates model output against actually-detected colonies in a subsequent window. The model is described formally in `ms/Point-process-model.Rmd` / `ms/Point-process-model_static.Rmd` — read these to understand the process/observation model before changing the fitting code. `j-foraging-distance.R` and `src/figures/queen-cells.R` are auxiliary analyses feeding this lineage.

Scripts suffixed `_pawsey` (e.g. `j-pp-static_pawsey.R`, `j-pp-static-get-data_pawsey.R`, `src/figures/j-pp-static-figures_pawsey.R`) are adapted for the Pawsey Setonix HPC cluster: they take a Slurm array index via `commandArgs()` to run one of several rolling 90-day time windows per job, rather than the single fixed window used interactively. Cluster job scripts live in `src/cluster/` (`a-setup-environment.sh` loads required modules — R, GDAL/GEOS/PROJ/UDUNITS, JAGS; `run_fit_JAGS.sh` / `run_fit_nimble.sh` are Slurm `sbatch` scripts for the `h-*` dynamic spatial occupancy fits).

JAGS/BUGS model definitions referenced by the `h-*` and `j-*` fitting scripts live as plain text in `src/model-files/`.

`src/figures/` mirrors the naming of the script that produces the data each figure script visualises (e.g. `f-multi-season-spatial-occupancy_figs.R` plots output from `f-multi-season-spatial-occupancy.R`).

## Key shared conventions

- `src/a-setup.R` defines `aggregate_data()` (spatial grid + temporal aggregation), `map_point_grid()` (leaflet visualisation of points + grid), and related helpers used across the `a`–`h` scripts. Read this file before modifying aggregation logic — it's depended on implicitly by filename convention, not by an explicit package API.
- `src/b-data-organisation.R` is the single source of truth for cleaning raw survey exports: parsing the surveillance `.xlsx`, separating colony records from point detections, deriving `presence`, `hive.removed`, distance-from-first-detection (`dist.0`), and other covariates. Any change here affects every downstream script.
- Distances for the point-process (`j-*`) lineage are computed in projected UTM 50S metres (EPSG:32750), not lat/long — reproject with `st_transform(crs = 32750)` before distance-based calculations, consistent with existing scripts.
- Fitted model objects and derived outputs (`.RData`, rasters, animations, figures) are written to `/out`, not version-controlled data — treat contents there as regenerable, not authoritative.
- Model reports/manuscripts live in `/ms` as R Markdown (`.Rmd`) with rendered `.html`/`.nb.html` outputs checked in alongside.
