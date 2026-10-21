# =============================================================================
# Parse Census of Governments: Finance (COG-F) Individual Unit Files
# Waves: 2012, 2017, 2022 — western states, local governments only
# =============================================================================
#
# Each wave has two fixed-width files (see the Technical Documentation PDF in
# each input/<year> COG-F/ folder):
#   - FinEstDAT: one row per government x item code (amount in $ thousands)
#   - Directory: one row per government (name, county, place, population)
#
# The record layout changed after 2012:
#
#                   2012 (Fin_GID)                2017+ (Fin_PID)
#   ID              14 chars                      12 chars
#   ID pos 1-2      COG-F internal state code     Census FIPS state code
#   County FIPS     directory pos 116-118         ID pos 4-6
#
# Both layouts are normalised here to Census FIPS (fips_state, fips_county),
# the key used for every downstream join (PAD-US, ACS, TIGER, RUCC).
# =============================================================================

library(tidyverse)
library(here)

cogf_files <- tribble(
  ~year, ~fin_file,                                        ~dir_file,
  2012,  "2012 COG-F/2012FinEstDAT_07242018modp_pu.txt",  "2012 COG-F/Fin_GID_2012.txt",
  2017,  "2017 COG-F/2017FinEstDAT_09202024modp_pu.txt",  "2017 COG-F/Fin_PID_2017.txt",
  2022,  "2022 COG-F/2022FinEstDAT_07152026modp.txt",     "2022 COG-F/Fin_PID_2022.txt"
)

western_state_fips <- c(
  "04" = "Arizona",    "06" = "California", "08" = "Colorado",
  "16" = "Idaho",      "30" = "Montana",    "32" = "Nevada",
  "35" = "New Mexico", "41" = "Oregon",     "49" = "Utah",
  "53" = "Washington", "56" = "Wyoming"
)

gov_types <- c(
  "1" = "County",
  "2" = "City",
  "3" = "Township",
  "4" = "Special District",
  "5" = "Independent School District or Educational Service Agency"
)

# =============================================================================
# RECORD LAYOUTS (1-based start/end positions from the technical documentation)
# =============================================================================

fin_layout <- list(
  "2012"  = fwf_cols(id = c(1, 14), item_code = c(15, 17),
                     amount_thousands = c(18, 29), flag = c(34, 34)),
  "2017+" = fwf_cols(id = c(1, 12), item_code = c(13, 15),
                     amount_thousands = c(16, 27), flag = c(32, 32))
)

dir_layout <- list(
  "2012"  = fwf_cols(id = c(1, 14), govt_name = c(15, 78), county_name = c(79, 113),
                     fips_state = c(114, 115), fips_county = c(116, 118),
                     fips_place = c(119, 123), population = c(124, 132)),
  "2017+" = fwf_cols(id = c(1, 12), govt_name = c(13, 76), county_name = c(77, 111),
                     fips_place = c(112, 116), population = c(117, 125))
)

# =============================================================================
# PARSE ONE WAVE
# =============================================================================

parse_cogf <- function(year, fin_file, dir_file) {
  layout <- if (year == 2012) "2012" else "2017+"
  message("Parsing COG-F ", year, " (", layout, " layout)...")

  fin <- read_fwf(here("input", fin_file), fin_layout[[layout]],
                  col_types = cols(amount_thousands = col_double(), .default = col_character()))

  dir <- read_fwf(here("input", dir_file), dir_layout[[layout]],
                  col_types = cols(population = col_double(), .default = col_character()))

  # 2017+ directories carry FIPS state/county only inside the ID
  if (layout == "2017+") {
    dir <- dir |> mutate(fips_state = substr(id, 1, 2), fips_county = substr(id, 4, 6))
  }

  fin |>
    inner_join(dir, by = "id") |>
    mutate(year = year, gov_type_code = substr(id, 3, 3)) |>
    filter(fips_state %in% names(western_state_fips),
           gov_type_code %in% names(gov_types))       # drops state governments
}

cogf_west <- pmap(cogf_files, parse_cogf) |>
  list_rbind() |>
  mutate(
    state_name    = western_state_fips[fips_state],
    gov_type_name = gov_types[gov_type_code],
    # Class codes by item code prefix (lookup: docs/Class Types.txt)
    class_type = case_when(
      substr(item_code, 1, 1) == "T"                                  ~ "R_Tax",
      substr(item_code, 1, 1) %in% c("B", "C", "D")                  ~ "R_IG",
      substr(item_code, 1, 1) == "A"                                  ~ "R_Ch",
      substr(item_code, 1, 1) == "U"                                  ~ "R_Oth",
      substr(item_code, 1, 1) %in% c("E", "I", "J", "F", "G", "K")  ~ "E_Dir",
      substr(item_code, 1, 1) %in% c("L", "M", "Q", "S")             ~ "E_IG",
      grepl("^[0-9]", item_code)                                      ~ "D"   # Debt
    )
  ) |>
  relocate(year, state_name, fips_state, fips_county, county_name, id, govt_name,
           gov_type_code, gov_type_name)

# --- Coverage check: governments and counties per state x wave ---------------
cogf_west |>
  group_by(year, state_name) |>
  summarise(n_govts    = n_distinct(id),
            n_counties = n_distinct(fips_county),
            .groups = "drop") |>
  pivot_wider(names_from = year, values_from = c(n_govts, n_counties)) |>
  print()

saveRDS(cogf_west, here("output", "cogf_west.rds"))
message("Saved: cogf_west.rds")
