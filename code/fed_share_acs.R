#=============================================================================
  # Spatial Analysis to Calculate Federal Share of County Lands in the West
  # Sources:
  #   - Protected Areas Database of the US (PAD-US)(shape file)
  #   - US Census TIGER/Line files (shape file)
  # Technical Documentation: PADUS_Standard_Tables_1.xlsx
  # =============================================================================# 

# Packages for reading, manipulating, and plotting shapefiles

library(sf)       # Read/write/manipulate vector spatial data (shapefiles, etc.)
library(terra)    # Raster and vector spatial data
library(ggplot2)  # Plotting, including geom_sf() for spatial data
library(here)

# --- Reading a shapefile ---
# shp <- sf::st_read("path/to/your/file.shp")

# --- Quick plot ---
# plot(shp)
# ggplot(shp) + geom_sf()

# --- Check coordinate reference system ---
# sf::st_crs(shp)

# --- Transform CRS (e.g., to WGS84) ---
# shp_wgs84 <- sf::st_transform(shp, crs = 4326)

library(sf)
library(dplyr)
library(tigris)
library(tidycensus)
options(tigris_use_cache = TRUE)

pad_dir <- here("input", "PAD-US") #Protected Areas of the US
gdb_path <- here("input", "SMA_WM.gdb")   #BLM Surface Management Agency data

# Western states: abbreviation, FIPS code, PAD-US file suffix
western_states <- tibble(
  abbrev = c("AZ", "CA", "CO", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY"),
  fips   = c("04", "06", "08", "16", "30", "32", "35", "41", "49", "53", "56")
)

# Function to calculate federal share of acres by county for one state
calc_fed_share <- function(abbrev, fips) {
  gdb_file <- file.path(pad_dir, paste0("PADUS4_1_State", abbrev, ".gdb"))
  comb_layer <- paste0("PADUS4_1Comb_DOD_Trib_NGP_Fee_Desig_Ease_State_", abbrev)

  if (!file.exists(gdb_file)) {
    message("Skipping ", abbrev, " — GDB not found: ", gdb_file)
    return(NULL)
  }

  message("Processing ", abbrev, "...")

  # Load combined PAD-US layer and filter to federal ownership
  pad_fed <- st_read(gdb_file, layer = comb_layer, quiet = TRUE) |>
    filter(Own_Type == "FED")

  # Load TIGER counties and reproject to match PAD-US CRS
  cty <- counties(state = fips, cb = TRUE, year = 2023) |>
    st_transform(st_crs(pad_fed))

  # Total county area from TIGER geometry
  cty_area <- cty |>
    mutate(total_area_m2 = as.numeric(st_area(geometry))) |>
    st_drop_geometry() |>
    select(GEOID, NAME, total_area_m2)

  # Intersect federal parcels with counties and sum area
  fed_area <- st_intersection(pad_fed, cty |> select(GEOID)) |>
    mutate(fed_area_m2 = as.numeric(st_area(Shape))) |>
    st_drop_geometry() |>
    group_by(GEOID) |>
    summarise(fed_area_m2 = sum(fed_area_m2), .groups = "drop")

  # Intersect USFS-managed parcels with counties and sum area
  usfs_area <- st_intersection(
    pad_fed |> filter(Mang_Name == "USFS"),
    cty |> select(GEOID)
  ) |>
    mutate(usfs_area_m2 = as.numeric(st_area(Shape))) |>
    st_drop_geometry() |>
    group_by(GEOID) |>
    summarise(usfs_area_m2 = sum(usfs_area_m2), .groups = "drop")

  # Intersect DOI-managed parcels (BLM, NPS, FWS, BOR) with counties and sum area
  doi_area <- st_intersection(
    pad_fed |> filter(Mang_Name %in% c("BLM", "NPS", "FWS", "BOR")),
    cty |> select(GEOID)
  ) |>
    mutate(doi_area_m2 = as.numeric(st_area(Shape))) |>
    st_drop_geometry() |>
    group_by(GEOID) |>
    summarise(doi_area_m2 = sum(doi_area_m2), .groups = "drop")

  # Join and compute share
  cty_area |>
    left_join(fed_area,  by = "GEOID") |>
    left_join(usfs_area, by = "GEOID") |>
    left_join(doi_area,  by = "GEOID") |>
    mutate(
      state             = abbrev,
      fed_area_m2       = coalesce(fed_area_m2,  0),
      usfs_area_m2      = coalesce(usfs_area_m2, 0),
      doi_area_m2       = coalesce(doi_area_m2,  0),
      fed_acres         = fed_area_m2        * 0.000247105,
      usfs_acres        = usfs_area_m2       * 0.000247105,
      doi_acres         = doi_area_m2        * 0.000247105,
      usfs_doi_acres    = (usfs_area_m2 + doi_area_m2) * 0.000247105,
      total_acres       = total_area_m2      * 0.000247105,
      # pmin(..., 100) caps at 100% — values slightly above 100 can occur when
      # PAD-US polygons extend marginally beyond TIGER county boundaries due to
      # geometric mismatches between the two source datasets.
      fed_share_pct      = pmin(round(fed_area_m2  / total_area_m2 * 100, 1), 100),
      usfs_share_pct     = pmin(round(usfs_area_m2 / total_area_m2 * 100, 1), 100),
      doi_share_pct      = pmin(round(doi_area_m2  / total_area_m2 * 100, 1), 100),
      usfs_doi_share_pct = pmin(round((usfs_area_m2 + doi_area_m2) / total_area_m2 * 100, 1), 100)
    ) |>
    select(state, GEOID, county = NAME, total_acres,
           fed_acres, usfs_acres, doi_acres, usfs_doi_acres,
           fed_share_pct, usfs_share_pct, doi_share_pct, usfs_doi_share_pct)
}

# Run for all western states and combine results
fed_share_west <- purrr::map2(
  western_states$abbrev,
  western_states$fips,
  calc_fed_share
) |>
  purrr::list_rbind()

# Preview results
fed_share_west |> arrange(state, desc(fed_share_pct))

# --- ACS county-level variables (2008–2012 5-year ACS) ---
# B19013_001: Median household income
# B11001_001: Total households
# B15003_001: Population 25+ (denominator for education attainment)
# B15003_022: Bachelor's degree
# B15003_023: Master's degree
# B15003_024: Professional school degree
# B15003_025: Doctorate degree
acs_raw <- get_acs(
  geography = "county",
  variables = c("B19013_001", "B11001_001",
                "B15003_001", "B15003_022", "B15003_023", "B15003_024", "B15003_025"),
  state     = western_states$fips,
  year      = 2012,
  survey    = "acs5",
  key       = Sys.getenv("CENSUS_API_KEY")
)

acs_income <- acs_raw |>
  filter(variable == "B19013_001") |>
  select(GEOID, med_hh_income = estimate, med_hh_income_moe = moe)

acs_households <- acs_raw |>
  filter(variable == "B11001_001") |>
  select(GEOID, n_households = estimate)

# Percent of population 25+ with a bachelor's degree or higher
acs_educ <- acs_raw |>
  filter(variable %in% c("B15003_001", "B15003_022", "B15003_023",
                         "B15003_024", "B15003_025")) |>
  select(GEOID, variable, estimate) |>
  pivot_wider(names_from = variable, values_from = estimate) |>
  mutate(
    pct_bachelors_plus = (B15003_022 + B15003_023 + B15003_024 + B15003_025) /
                         B15003_001 * 100
  ) |>
  select(GEOID, pct_bachelors_plus)

# Join income, household count, and education onto federal land share data
fed_share_west <- fed_share_west |>
  left_join(acs_income,     by = "GEOID") |>
  left_join(acs_households, by = "GEOID") |>
  left_join(acs_educ,       by = "GEOID")

# Save results
write.csv(fed_share_west, here("output", "fed_share_west.csv"), row.names = FALSE)
message("Saved: fed_share_west.csv")

# Join federal land share onto COG-F microdata by state + county FIPS
cogf_2012w <- readRDS(here("output", "cogf_2012w.rds"))

cogf_2012w_fed <- cogf_2012w |>
  left_join(
    fed_share_west |>
      mutate(
        fips_state  = substr(GEOID, 1, 2),
        fips_county = substr(GEOID, 3, 5)
      ) |>
      select(fips_state, fips_county, fed_share_pct, usfs_share_pct, doi_share_pct, usfs_doi_share_pct, med_hh_income, med_hh_income_moe),
    by = c("fips_state", "fips_county")
  )

saveRDS(cogf_2012w_fed, here("output", "cogf_2012w_fed.rds"))
message("Saved: cogf_2012w_fed.rds")

write.xlsx(cogf_2012w_fed, here("output", "cogf_2012w_fed.xlsx"), rowNames = FALSE)
