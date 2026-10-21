
library(tidyverse)
library(openxlsx)
library(readxl)
library(here)
# tigris and sf are no longer needed here — spatial computation has been moved
# to unincorp_share.R, which writes unincorp_share.csv read below.
#-----------------------------------------------------------------------------------------------------
# Note that rowSums(..., na.rm = TRUE) treats any NA class codes as zero in the totals — worth keeping in mind if a county is missing a category entirely.
#-----------------------------------------------------------------------------------------------------

data1 <- readRDS(here("R Code and Data", "cogf_2012w_fed.rds"))

# 1. Count of distinct governments per county
gov_counts <- data1 |>
  distinct(state_code, county_name, id) |>
  count(state_code, county_name, name = "n_govts")

# 2. Aggregate amount_thousands by class_type (Class Code) and county
cty_class_agg <- data1 |>
  filter(!is.na(class_type)) |>
  group_by(state_code, county_name, class_type) |>
  summarise(total_thousands = sum(amount_thousands, na.rm = TRUE), .groups = "drop")

# Federal IG Revenue only (item codes starting with "B")
fed_ig_agg <- data1 |>
  filter(substr(item_code, 1, 1) == "B") |>
  group_by(state_code, county_name) |>
  summarise(R_IG_Fed = sum(amount_thousands, na.rm = TRUE), .groups = "drop")

# 3. Join government counts and state_name onto aggregated data
state_names <- data1 |>
  distinct(state_code, state_name)

# Pull year, FIPS codes, and population from county government rows (gov_type_code == "1").
# For consolidated city-counties (e.g. San Francisco, Denver), no county government row
# exists; fall back to the city government row (gov_type_code == "2").
cty_meta <- data1 |>
  filter(gov_type_code %in% c("1", "2")) |>
  distinct(state_code, county_name, gov_type_code, year, fips_county, fips_place, population) |>
  arrange(state_code, county_name, gov_type_code) |>   # "1" sorts before "2"
  group_by(state_code, county_name) |>
  slice(1) |>                                          # keep county govt; city only if no county
  ungroup() |>
  select(-gov_type_code)


data2 <- cty_class_agg |>
  left_join(gov_counts,  by = c("state_code", "county_name")) |>
  left_join(state_names, by = "state_code") |>
  left_join(cty_meta,    by = c("state_code", "county_name")) |>
  pivot_wider(names_from = class_type, values_from = total_thousands) |>
  left_join(fed_ig_agg,  by = c("state_code", "county_name")) |>
  select(state_code, state_name, county_name, year, fips_county, fips_place, population,
         n_govts, R_Tax, R_IG, R_IG_Fed, R_Ch, R_Oth, E_Dir, E_IG, D) |>
  mutate(
    Total_Revenue        = rowSums(across(c(R_Tax, R_IG, R_Ch, R_Oth)), na.rm = TRUE),
    Total_Expenditures   = rowSums(across(c(E_Dir, E_IG)), na.rm = TRUE),
    PerCap_TaxBurden     = (R_Tax + R_Ch) / population,
    PerCap_Expenditures  = E_Dir / population,
    GovtsPer10k          = n_govts / population * 10000
  )

# Join federal land share (%) from shapefile_work.R via GEOID (state FIPS + county FIPS)
state_fips_lookup <- tibble(
  state_name = c("Arizona", "California", "Colorado", "Idaho",    "Montana",
                 "Nevada",  "New Mexico", "Oregon",   "Utah",     "Washington", "Wyoming"),
  state_fips = c("04",      "06",         "08",       "16",       "30",
                 "32",      "35",         "41",       "49",       "53",         "56")
)

fed_share_join <- read.csv(here("R Code and Data", "fed_share_west.csv"),
                           colClasses = c(GEOID = "character")) |>
  select(GEOID, fed_share_pct, usfs_share_pct, doi_share_pct, usfs_doi_share_pct,
         med_hh_income, med_hh_income_moe, n_households, pct_bachelors_plus)

# Download ERS Rural-Urban Continuum Codes 2013 (metro = codes 1-3, non-metro = codes 4-9)
rucc_path <- here("R Code and Data", "ruralurbancodes2013.xls")
if (!file.exists(rucc_path)) {
  download.file(
    "https://www.ers.usda.gov/media/5769/2013-rural-urban-continuum-codes.xls",
    destfile = rucc_path, mode = "wb"
  )
}

metro_join <- read_xls(rucc_path) |>
  select(GEOID = FIPS, RUCC_2013) |>
  mutate(
    GEOID = as.character(GEOID),
    metro = if_else(RUCC_2013 <= 3, "Metro", "Non-Metro")
  ) |>
  select(GEOID, RUCC_2013, metro)

# =============================================================================
# Read pre-computed unincorporated area share (one-time spatial computation).
# Run unincorp_share.R to regenerate this file if county or place boundaries
# need to be updated. GEOID is the 5-digit Census FIPS (state + county).
# =============================================================================

unincorp_path <- here("R Code and Data", "unincorp_share.csv")
if (!file.exists(unincorp_path)) {
  stop("unincorp_share.csv not found. Run unincorp_share.R first to generate it.")
}
unincorp_join <- read_csv(unincorp_path,
                          col_types = cols(GEOID = col_character()),
                          show_col_types = FALSE)

# =============================================================================

data2 <- data2 |>
  left_join(state_fips_lookup, by = "state_name") |>
  mutate(GEOID = if_else(!is.na(state_fips), paste0(state_fips, fips_county), NA_character_)) |>
  left_join(fed_share_join,  by = "GEOID") |>
  left_join(metro_join,      by = "GEOID") |>
  left_join(unincorp_join,   by = "GEOID") |>
  select(-state_fips, -GEOID) |>
  # HH burden metrics: per-household tax/expenditure amounts as a share of
  # median household income. Steps: (1) convert from $thousands to dollars
  # (* 1000), (2) divide by household count for per-HH amount, (3) divide by
  # median HH income and multiply by 100 to express as a percentage.
  mutate(
    HH_TaxBurden    = (R_Tax + R_Ch) * 1000 / n_households, 
    HH_Expenditures = E_Dir          * 1000 / n_households 
  ) |>
  relocate(c(usfs_share_pct, doi_share_pct, usfs_doi_share_pct, unincorp_share_pct, med_hh_income, med_hh_income_moe), .after = fed_share_pct) |>
  relocate(c(metro, RUCC_2013), .after = n_govts) |>
  relocate(c(HH_TaxBurden, HH_Expenditures), .after = PerCap_Expenditures)

# Save as .rds
saveRDS(data2, here("R Code and Data", "Cty_Agg.rds"))
message("Saved: Cty_Agg.rds")

# Also save as CSV if a plain-text copy is needed
write.xlsx(data2, here("R Code and Data", "County_Aggregate.xlsx"))

# Save as xlsx with full column names for readability
data2_xlsx <- data2 |>
  rename(
    "Number of Governments"   = n_govts,
    "Govts per 10,000 Pop"    = GovtsPer10k,
    "Number of Households"    = n_households,
    "General Revenue-Taxes"   = R_Tax,
    "General Revenue-IGR"     = R_IG,
    "General Revenue-IGR Fed" = R_IG_Fed,
    "General Revenue-Charges" = R_Ch,
    "General Revenue-Misc"    = R_Oth,
    "Expenditures-Direct"     = E_Dir,
    "Expenditures-IG"         = E_IG,
    "Debt"                    = D,
    "Total Revenue"           = Total_Revenue,
    "Total Expenditures"      = Total_Expenditures,
    "Per Capita Tax Burden"   = PerCap_TaxBurden,
    "Per Capita Expenditures" = PerCap_Expenditures,
    "Metro/Non-Metro"         = metro,
    "RUCC Code (2013)"        = RUCC_2013,
    "Fed Share (%)"           = fed_share_pct,
    "USFS Share (%)"          = usfs_share_pct,
    "DOI Share (%)"           = doi_share_pct,
    "USFS DOI Share (%)"      = usfs_doi_share_pct,
    "Median Household Income" = med_hh_income,
    "Median HH Income MOE"    = med_hh_income_moe,
    "Bachelor's Degree or Higher (%)" = pct_bachelors_plus,
    "Unincorporated Share (%)"        = unincorp_share_pct,
    "HH Tax Burden (% of Med HH Income)"   = HH_TaxBurden,
    "HH Expenditures (% of Med HH Income)" = HH_Expenditures,
  )

write.xlsx(data2_xlsx, file = here("R Code and Data", "County_Aggregate.xlsx"),
           sheetName = "Cty_Agg", rowNames = FALSE)
message("Saved: County_Aggregate.xlsx")

