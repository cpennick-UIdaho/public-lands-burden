# =============================================================================
# County-Level Aggregation — 2022 COG-F
# =============================================================================
# Reads cogf_2022w.rds, aggregates to county level, pulls 2017 ACS variables,
# joins federal land share (fed_land_share.csv), RUCC metro codes, and
# unincorporated area share, then saves Cty_Agg_2022.rds.
#
# Differences from Cty_Agg.R (2012):
#   - Input: cogf_2022w.rds (from All_Gov_West_2022.R)
#   - ACS year: 2022 5-year ACS → acs_cache_2022.csv
#   - GEOID: built directly from fips_state + fips_county — 2017 data uses
#     Census FIPS in the ID; no COG-F internal code lookup needed
#   - D (debt) class may be absent — handled with any_of() in select()
#
# One-time inputs (must exist before running):
#   fed_land_share.csv   — run fed_land_share.R once to generate
#   unincorp_share.csv   — run unincorp_share.R once to generate
#
# ACS pull is gated: if acs_cache_2022.csv exists it is loaded from disk.
# Delete it to force a fresh pull from the Census API.
# =============================================================================

library(tidyverse)
library(openxlsx)
library(readxl)
library(tidycensus)

# Explicit repo root — all file paths resolve here regardless of active .Rproj
repo <- "C:/Users/chelseap/OneDrive - University of Idaho/Documents/GitHub/public-lands-burden"

# =============================================================================
# PART 1 — COG-F COUNTY AGGREGATION
# =============================================================================

data1 <- readRDS(file.path(repo, "cogf_2022w.rds"))

gov_counts <- data1 |>
  distinct(state_code, county_name, id) |>
  count(state_code, county_name, name = "n_govts")

cty_class_agg <- data1 |>
  filter(!is.na(class_type)) |>
  group_by(state_code, county_name, class_type) |>
  summarise(total_thousands = sum(amount_thousands, na.rm = TRUE), .groups = "drop")

fed_ig_agg <- data1 |>
  filter(substr(item_code, 1, 1) == "B") |>
  group_by(state_code, county_name) |>
  summarise(R_IG_Fed = sum(amount_thousands, na.rm = TRUE), .groups = "drop")

state_names <- data1 |>
  distinct(state_code, state_name)

# Pull year, FIPS codes, and population from county government rows.
# Fall back to city row for consolidated city-counties.
cty_meta <- data1 |>
  filter(gov_type_code %in% c("1", "2")) |>
  distinct(state_code, county_name, gov_type_code, year,
           fips_state, fips_county, fips_place, population) |>
  arrange(state_code, county_name, gov_type_code) |>
  group_by(state_code, county_name) |>
  slice(1) |>
  ungroup() |>
  select(-gov_type_code)

data2 <- cty_class_agg |>
  left_join(gov_counts,  by = c("state_code", "county_name")) |>
  left_join(state_names, by = "state_code") |>
  left_join(cty_meta,    by = c("state_code", "county_name")) |>
  pivot_wider(names_from = class_type, values_from = total_thousands) |>
  left_join(fed_ig_agg,  by = c("state_code", "county_name")) |>
  # Use any_of() — D (debt) may be absent in some COG-F years
  select(state_code, state_name, county_name, year, fips_state, fips_county,
         fips_place, population, n_govts,
         any_of(c("R_Tax", "R_IG", "R_IG_Fed", "R_Ch", "R_Oth", "E_Dir", "E_IG", "D"))) |>
  mutate(
    Total_Revenue       = rowSums(across(any_of(c("R_Tax", "R_IG", "R_Ch", "R_Oth"))),
                                  na.rm = TRUE),
    Total_Expenditures  = rowSums(across(any_of(c("E_Dir", "E_IG"))), na.rm = TRUE),
    PerCap_TaxBurden    = (R_Tax + R_Ch) / population,
    PerCap_Expenditures = E_Dir / population,
    GovtsPer10k         = n_govts / population * 10000
  )

# =============================================================================
# PART 2 — BUILD GEOID (2017: fips_state is already Census FIPS; no lookup)
# =============================================================================

data2 <- data2 |>
  mutate(GEOID = paste0(fips_state, fips_county))

# =============================================================================
# PART 3 — ACS PULL (gated; 2022 5-year ACS)
# =============================================================================

acs_cache_path <- file.path(repo, "acs_cache_2022.csv")
western_fips   <- unique(data1$fips_state)

if (file.exists(acs_cache_path)) {
  message("acs_cache_2022.csv found — loading from disk, skipping ACS API call.")
  acs_join <- read_csv(acs_cache_path,
                       col_types = cols(GEOID = col_character()),
                       show_col_types = FALSE)
} else {
  message("Pulling 2022 ACS 5-year estimates...")

  acs_raw <- get_acs(
    geography = "county",
    variables = c(
      "B19013_001",
      "B11001_001",
      "B15003_001",
      "B15003_022", "B15003_023", "B15003_024", "B15003_025"
    ),
    state  = western_fips,
    year   = 2022,
    survey = "acs5",
    key    = Sys.getenv("CENSUS_API_KEY")
  )

  acs_income <- acs_raw |>
    filter(variable == "B19013_001") |>
    select(GEOID, med_hh_income = estimate, med_hh_income_moe = moe)

  acs_households <- acs_raw |>
    filter(variable == "B11001_001") |>
    select(GEOID, n_households = estimate)

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

  acs_join <- acs_income |>
    left_join(acs_households, by = "GEOID") |>
    left_join(acs_educ,       by = "GEOID")

  write_csv(acs_join, acs_cache_path)
  message("Saved: acs_cache_2022.csv — ", nrow(acs_join), " counties")
}

# =============================================================================
# PART 4 — JOIN SUPPLEMENTARY DATA
# =============================================================================

land_share_path <- file.path(repo, "fed_land_share.csv")
if (!file.exists(land_share_path)) {
  stop("fed_land_share.csv not found. Run fed_land_share.R first to generate it.")
}
fed_land_share <- read.csv(land_share_path, colClasses = c(GEOID = "character")) |>
  select(GEOID, fed_share_pct, usfs_share_pct, doi_share_pct, usfs_doi_share_pct)

rucc_path <- file.path(repo, "ruralurbancodes2013.xls")
if (!file.exists(rucc_path)) {
  download.file("https://www.ers.usda.gov/media/5769/2013-rural-urban-continuum-codes.xls",
                destfile = rucc_path, mode = "wb")
}
metro_join <- read_xls(rucc_path) |>
  select(GEOID = FIPS, RUCC_2013) |>
  mutate(GEOID = as.character(GEOID),
         metro = if_else(RUCC_2013 <= 3, "Metro", "Non-Metro")) |>
  select(GEOID, RUCC_2013, metro)

unincorp_path <- file.path(repo, "unincorp_share.csv")
if (!file.exists(unincorp_path)) {
  stop("unincorp_share.csv not found. Run unincorp_share.R first to generate it.")
}
unincorp_join <- read_csv(unincorp_path,
                          col_types = cols(GEOID = col_character()),
                          show_col_types = FALSE)

data2 <- data2 |>
  left_join(fed_land_share, by = "GEOID") |>
  left_join(acs_join,       by = "GEOID") |>
  left_join(metro_join,     by = "GEOID") |>
  left_join(unincorp_join,  by = "GEOID") |>
  select(-GEOID) |>
  mutate(
    HH_TaxBurden    = (R_Tax + R_Ch) * 1000 / n_households,
    HH_Expenditures = E_Dir          * 1000 / n_households
  ) |>
  relocate(c(usfs_share_pct, doi_share_pct, usfs_doi_share_pct, unincorp_share_pct,
             med_hh_income, med_hh_income_moe), .after = fed_share_pct) |>
  relocate(c(metro, RUCC_2013), .after = n_govts) |>
  relocate(c(HH_TaxBurden, HH_Expenditures), .after = PerCap_Expenditures)

# =============================================================================
# PART 5 — SAVE
# =============================================================================

saveRDS(data2, file.path(repo, "Cty_Agg_2022.rds"))
message("Saved: Cty_Agg_2022.rds — ", nrow(data2), " counties")

data2_xlsx <- data2 |>
  rename(
    "Number of Governments"             = n_govts,
    "Govts per 10,000 Pop"              = GovtsPer10k,
    "Number of Households"              = n_households,
    "General Revenue-Taxes"             = R_Tax,
    "General Revenue-IGR"               = R_IG,
    "General Revenue-IGR Fed"           = R_IG_Fed,
    "General Revenue-Charges"           = R_Ch,
    "General Revenue-Misc"              = R_Oth,
    "Expenditures-Direct"               = E_Dir,
    "Expenditures-IG"                   = E_IG,
    "Total Revenue"                     = Total_Revenue,
    "Total Expenditures"                = Total_Expenditures,
    "Per Capita Tax Burden"             = PerCap_TaxBurden,
    "Per Capita Expenditures"           = PerCap_Expenditures,
    "Metro/Non-Metro"                   = metro,
    "RUCC Code (2013)"                  = RUCC_2013,
    "Fed Share (%)"                     = fed_share_pct,
    "USFS Share (%)"                    = usfs_share_pct,
    "DOI Share (%)"                     = doi_share_pct,
    "USFS DOI Share (%)"                = usfs_doi_share_pct,
    "Median Household Income"           = med_hh_income,
    "Median HH Income MOE"              = med_hh_income_moe,
    "Bachelor's Degree or Higher (%)"   = pct_bachelors_plus,
    "Unincorporated Share (%)"          = unincorp_share_pct,
    "HH Tax Burden ($)"                 = HH_TaxBurden,
    "HH Expenditures ($)"               = HH_Expenditures
  )

write.xlsx(data2_xlsx,
           file      = file.path(repo,"County_Aggregate_2022.xlsx"),
           sheetName = "Cty_Agg_2022",
           rowNames  = FALSE)
message("Saved: County_Aggregate_2022.xlsx")
