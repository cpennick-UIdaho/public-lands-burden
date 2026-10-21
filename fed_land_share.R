# =============================================================================
# PAD-US Spatial Intersection — Federal Land Share by County
# =============================================================================
# ONE-TIME computation. Re-run only if the PAD-US version changes.
# Output: fed_land_share.csv  (one row per county; purely spatial, no ACS)
#
# To force a recompute: delete fed_land_share.csv and source this script.
# All downstream scripts (fed_share_acs.R, Cty_Agg*.R) read the CSV; they
# never need to touch PAD-US GDB files directly.
#
# Sources:
#   PAD-US v4.1 state GDB files in DATA/PAD-US/
#   TIGER/Line county cartographic boundaries (year = 2023) via tigris
# =============================================================================

library(sf)
library(dplyr)
library(purrr)
library(tigris)
library(here)

options(tigris_use_cache = TRUE)

pad_dir <- here("DATA", "PAD-US")

western_states <- tibble(
  abbrev = c("AZ", "CA", "CO", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY"),
  fips   = c("04", "06", "08", "16", "30", "32", "35", "41", "49", "53", "56")
)

# =============================================================================
# Per-state intersection function
# =============================================================================
# Returns one row per county with total, USFS, DOI, and combined land share (%).
# pmin(..., 100) caps at 100% — slight overshoots occur when PAD-US polygons
# extend marginally beyond TIGER county boundaries (geometric mismatch).

calc_fed_share <- function(abbrev, fips) {

  gdb_file   <- file.path(pad_dir, paste0("PADUS4_1_State", abbrev, ".gdb"))
  comb_layer <- paste0("PADUS4_1Comb_DOD_Trib_NGP_Fee_Desig_Ease_State_", abbrev)

  if (!file.exists(gdb_file)) {
    message("Skipping ", abbrev, " — GDB not found: ", gdb_file)
    return(NULL)
  }

  message("Processing ", abbrev, "...")

  pad_fed <- st_read(gdb_file, layer = comb_layer, quiet = TRUE) |>
    filter(Own_Type == "FED")

  cty <- counties(state = fips, cb = TRUE, year = 2023) |>
    st_transform(st_crs(pad_fed))

  cty_area <- cty |>
    mutate(total_area_m2 = as.numeric(st_area(geometry))) |>
    st_drop_geometry() |>
    select(GEOID, NAME, total_area_m2)

  # All federal land
  fed_area <- st_intersection(pad_fed, cty |> select(GEOID)) |>
    mutate(fed_area_m2 = as.numeric(st_area(Shape))) |>
    st_drop_geometry() |>
    group_by(GEOID) |>
    summarise(fed_area_m2 = sum(fed_area_m2), .groups = "drop")

  # USFS only
  usfs_area <- st_intersection(
    pad_fed |> filter(Mang_Name == "USFS"),
    cty |> select(GEOID)
  ) |>
    mutate(usfs_area_m2 = as.numeric(st_area(Shape))) |>
    st_drop_geometry() |>
    group_by(GEOID) |>
    summarise(usfs_area_m2 = sum(usfs_area_m2), .groups = "drop")

  # DOI agencies (BLM, NPS, FWS, BOR)
  doi_area <- st_intersection(
    pad_fed |> filter(Mang_Name %in% c("BLM", "NPS", "FWS", "BOR")),
    cty |> select(GEOID)
  ) |>
    mutate(doi_area_m2 = as.numeric(st_area(Shape))) |>
    st_drop_geometry() |>
    group_by(GEOID) |>
    summarise(doi_area_m2 = sum(doi_area_m2), .groups = "drop")

  cty_area |>
    left_join(fed_area,  by = "GEOID") |>
    left_join(usfs_area, by = "GEOID") |>
    left_join(doi_area,  by = "GEOID") |>
    mutate(
      state          = abbrev,
      fed_area_m2    = coalesce(fed_area_m2,  0),
      usfs_area_m2   = coalesce(usfs_area_m2, 0),
      doi_area_m2    = coalesce(doi_area_m2,  0),
      fed_acres         = fed_area_m2                      * 0.000247105,
      usfs_acres        = usfs_area_m2                     * 0.000247105,
      doi_acres         = doi_area_m2                      * 0.000247105,
      usfs_doi_acres    = (usfs_area_m2 + doi_area_m2)     * 0.000247105,
      total_acres       = total_area_m2                    * 0.000247105,
      fed_share_pct      = pmin(round(fed_area_m2  / total_area_m2 * 100, 1), 100),
      usfs_share_pct     = pmin(round(usfs_area_m2 / total_area_m2 * 100, 1), 100),
      doi_share_pct      = pmin(round(doi_area_m2  / total_area_m2 * 100, 1), 100),
      usfs_doi_share_pct = pmin(round(
        (usfs_area_m2 + doi_area_m2) / total_area_m2 * 100, 1), 100)
    ) |>
    select(state, GEOID, county = NAME, total_acres,
           fed_acres, usfs_acres, doi_acres, usfs_doi_acres,
           fed_share_pct, usfs_share_pct, doi_share_pct, usfs_doi_share_pct)
}

# =============================================================================
# Run (gated — skipped if output already exists)
# =============================================================================

out_path <- here("fed_land_share.csv")

if (file.exists(out_path)) {
  message("fed_land_share.csv already exists — skipping PAD-US intersection.")
  message("Delete the file and re-run this script to recompute.")
} else {
  message("Running PAD-US spatial intersection for 11 western states...")
  message("This typically takes 10–20 minutes depending on hardware.")

  fed_land_share <- map2(western_states$abbrev, western_states$fips,
                         calc_fed_share) |>
    list_rbind()

  write.csv(fed_land_share, out_path, row.names = FALSE)
  message("Saved: fed_land_share.csv — ", nrow(fed_land_share), " counties")
}
