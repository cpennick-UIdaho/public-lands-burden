
# =============================================================================
# unincorp_share.R
# =============================================================================
# Computes unincorporated area share (%) for each western US county.
#
# Method: intersect 2013 TIGER county boundaries with 2020 TIGER incorporated
# place boundaries (CLASSFP C1/C2/C6). For each county, unincorporated share =
# (county area - incorporated place area) / county area * 100.
#
# Output: unincorp_share.csv  (GEOID, unincorp_share_pct)
#         One row per county; GEOID is the 5-digit Census FIPS (state + county).
#
# Run this script once and commit the CSV. Cty_Agg.R reads the CSV directly,
# avoiding the expensive spatial download and intersection on every run.
#
# Notes on vintage mismatches:
#   - County boundaries: 2013 cartographic boundary files (earliest available;
#     effectively identical to 2012).
#   - Incorporated place boundaries: 2020 cb files (2013 state-level place files
#     are not hosted by Census; place boundaries change slowly — acceptable
#     substitution for a 2012 cross-section).
# =============================================================================

library(tidyverse)
library(sf)
library(tigris)
library(here)

options(tigris_use_cache = TRUE)

# Western states Census FIPS codes
western_state_fips <- c("04", "06", "08", "16", "30", "32", "35", "41",
                         "49", "53", "56")

# --- Download TIGER county boundaries (2013) ----------------------------------
message("Downloading TIGER county boundaries (2013)...")
counties_tiger <- map(western_state_fips,
                      \(s) counties(state = s, year = 2013, cb = TRUE)) |>
  bind_rows() |>
  st_transform(5070) |>       # NAD83 / Conus Albers (equal-area projection)
  mutate(county_area = as.numeric(st_area(geometry)))

message("  ", nrow(counties_tiger), " counties loaded.")

# --- Download TIGER incorporated places (2020) --------------------------------
message("Downloading TIGER incorporated places (2020)...")
places_tiger <- map(western_state_fips,
                    \(s) places(state = s, year = 2020, cb = FALSE)) |>
  bind_rows() |>
  st_transform(5070) |>
  filter(CLASSFP %in% c("C1", "C2", "C6"))   # Active incorporated municipalities

message("  ", nrow(places_tiger), " incorporated places loaded.")

# --- Intersect places with counties -------------------------------------------
message("Computing place-county intersections (may take several minutes)...")
place_cty_int <- st_intersection(
  places_tiger  |> select(geometry),
  counties_tiger |> select(GEOID, county_area)
)

# --- Summarise incorporated area per county -----------------------------------
incorp_by_cty <- place_cty_int |>
  mutate(incorp_area = as.numeric(st_area(geometry))) |>
  st_drop_geometry() |>
  group_by(GEOID) |>
  summarise(total_incorp_area = sum(incorp_area), .groups = "drop")

# --- Compute unincorporated share ---------------------------------------------
unincorp_share <- counties_tiger |>
  st_drop_geometry() |>
  select(GEOID, county_area) |>
  left_join(incorp_by_cty, by = "GEOID") |>
  mutate(
    total_incorp_area  = replace_na(total_incorp_area, 0),
    unincorp_share_pct = (county_area - total_incorp_area) / county_area * 100
  ) |>
  select(GEOID, unincorp_share_pct)

message("Unincorporated area shares calculated for ", nrow(unincorp_share),
        " counties.")

# --- Save ---------------------------------------------------------------------
out_path <- here("output", "unincorp_share.csv")
write_csv(unincorp_share, out_path)
message("Saved: ", out_path)
