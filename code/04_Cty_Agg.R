
library(tidyverse)
library(openxlsx)
library(readxl)
library(here)
# tigris and sf are no longer needed here — spatial computation has been moved
# to 03_unincorp_share.R, which writes unincorp_share.csv read below.
#-----------------------------------------------------------------------------------------------------
# Note that rowSums(..., na.rm = TRUE) treats any NA class codes as zero in the totals — worth keeping in mind if a county is missing a category entirely.
#-----------------------------------------------------------------------------------------------------

data1 <- readRDS(here("output", "cogf_west.rds")) |>
  mutate(GEOID = paste0(fips_state, fips_county))

# One row per county x COG-F wave
cty_keys <- c("year", "GEOID")

# 1. Count of distinct governments per county
gov_counts <- data1 |>
  distinct(year, GEOID, id) |>
  count(year, GEOID, name = "n_govts")

# 2. Aggregate amount_thousands by class_type (Class Code) and county
cty_class_agg <- data1 |>
  filter(!is.na(class_type)) |>
  group_by(year, GEOID, class_type) |>
  summarise(total_thousands = sum(amount_thousands, na.rm = TRUE), .groups = "drop")

# Federal IG Revenue only (item codes starting with "B")
fed_ig_agg <- data1 |>
  filter(substr(item_code, 1, 1) == "B") |>
  group_by(year, GEOID) |>
  summarise(R_IG_Fed = sum(amount_thousands, na.rm = TRUE), .groups = "drop")

# Pull names, FIPS place, and population from county government rows (gov_type_code == "1").
# For consolidated city-counties (e.g. San Francisco, Denver), no county government row
# exists; fall back to the city government row (gov_type_code == "2").
cty_meta <- data1 |>
  filter(gov_type_code %in% c("1", "2")) |>
  distinct(year, GEOID, gov_type_code, state_name, fips_county, county_name, fips_place, population) |>
  arrange(year, GEOID, gov_type_code) |>   # "1" sorts before "2"
  group_by(year, GEOID) |>
  slice(1) |>                              # keep county govt; city only if no county
  ungroup() |>
  select(-gov_type_code)


data2 <- cty_class_agg |>
  pivot_wider(names_from = class_type, values_from = total_thousands) |>
  left_join(gov_counts, by = cty_keys) |>
  left_join(cty_meta,   by = cty_keys) |>
  left_join(fed_ig_agg, by = cty_keys) |>
  select(year, GEOID, state_name, county_name, fips_county, fips_place, population,
         n_govts, R_Tax, R_IG, R_IG_Fed, R_Ch, R_Oth, E_Dir, E_IG, D) |>
  mutate(
    Total_Revenue        = rowSums(across(c(R_Tax, R_IG, R_Ch, R_Oth)), na.rm = TRUE),
    Total_Expenditures   = rowSums(across(c(E_Dir, E_IG)), na.rm = TRUE),
    PerCap_TaxBurden     = (R_Tax + R_Ch) / population,
    PerCap_Expenditures  = E_Dir / population,
    GovtsPer10k          = n_govts / population * 10000
  )

# Federal land share (%) from 02_fed_share_acs.R — single PAD-US vintage, joined to every wave
fed_share_join <- read.csv(here("output", "fed_share_west.csv"),
                           colClasses = c(GEOID = "character")) |>
  select(GEOID, fed_share_pct, usfs_share_pct, doi_share_pct, usfs_doi_share_pct)

# ACS 5-year income, households, and education from 02_fed_share_acs.R — matched to each wave
acs_join <- read.csv(here("output", "acs_west.csv"),
                     colClasses = c(GEOID = "character"))

# ERS Rural-Urban Continuum Codes 2013 (metro = codes 1-3, non-metro = codes 4-9).
# Held fixed across waves so metro status is constant within county.
rucc_path <- here("input", "ruralurbancodes2013.xls")
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
# Run 03_unincorp_share.R to regenerate this file if county or place boundaries
# need to be updated. GEOID is the 5-digit Census FIPS (state + county).
# =============================================================================

unincorp_path <- here("output", "unincorp_share.csv")
if (!file.exists(unincorp_path)) {
  stop("unincorp_share.csv not found. Run 03_unincorp_share.R first to generate it.")
}
unincorp_join <- read_csv(unincorp_path,
                          col_types = cols(GEOID = col_character()),
                          show_col_types = FALSE)

# =============================================================================

data2 <- data2 |>
  left_join(fed_share_join, by = "GEOID") |>
  left_join(acs_join,       by = cty_keys) |>
  left_join(metro_join,     by = "GEOID") |>
  left_join(unincorp_join,  by = "GEOID") |>
  # HH burden metrics: nominal dollars per household. Steps: (1) convert from
  # $thousands to dollars (* 1000), (2) divide by ACS household count.
  mutate(
    HH_TaxBurden    = (R_Tax + R_Ch) * 1000 / n_households,
    HH_Expenditures = E_Dir          * 1000 / n_households
  ) |>
  relocate(c(usfs_share_pct, doi_share_pct, usfs_doi_share_pct, unincorp_share_pct, med_hh_income, med_hh_income_moe), .after = fed_share_pct) |>
  relocate(c(metro, RUCC_2013), .after = n_govts) |>
  relocate(c(HH_TaxBurden, HH_Expenditures), .after = PerCap_Expenditures) |>
  arrange(GEOID, year)

stopifnot(!anyDuplicated(data2[cty_keys]))

# Save as .rds
saveRDS(data2, here("output", "Cty_Agg.rds"))
message("Saved: Cty_Agg.rds")

# Also save as CSV if a plain-text copy is needed
write.xlsx(data2, here("output", "County_Aggregate.xlsx"))

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
    "HH Tax Burden ($ per Household)"     = HH_TaxBurden,
    "HH Expenditures ($ per Household)"   = HH_Expenditures,
  )

write.xlsx(data2_xlsx, file = here("output", "County_Aggregate.xlsx"),
           sheetName = "Cty_Agg", rowNames = FALSE)
message("Saved: County_Aggregate.xlsx")

