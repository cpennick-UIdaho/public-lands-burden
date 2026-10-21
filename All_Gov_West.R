# =============================================================================
# Parse 2012 State & Local Government Finance Data
# Sources:
#   - 2012FinEstDAT_07242018modp_pu.txt  (individual unit finance data)
#   - Fin_GID_2012.txt                   (government identifier / directory)
# Technical Documentation: 2012_S_L_Indiv_Unit_Data_File_Tech_Doc.pdf
# =============================================================================
#
# FINANCE FILE record layout (34 chars, fixed-width, no delimiter):
#   Pos  1-14 : ID (14 chars)
#     Pos  1-2 : State code
#     Pos  3   : Type-of-government code
#     Pos  4-6 : County FIPS
#     Pos  7-9 : Unit identifier
#     Pos 10-14: Should be "00000" (unit not part of another govt)
#   Pos 15-17  : Item code (3 chars)
#   Pos 18-29  : Amount in thousands of dollars (12 chars, numeric)
#   Pos 30-33  : Year of data (4 chars)
#   Pos 34     : Imputation type / item data flag (1 char): "I"=Imputed, "R"=Reported
#
# GID FILE record layout (153 chars, fixed-width, Windows line endings \r\n):
#   Pos   1-14  : ID code (14 chars) — join key to finance file
#   Pos  15-78  : Government name (64 chars)
#   Pos  79-113 : County name (35 chars)
#   Pos 114-115 : FIPS state code (2 chars)
#   Pos 116-118 : FIPS county code (3 chars)
#   Pos 119-123 : FIPS place code (5 chars)
#   Pos 124-132 : Population (9 chars, numeric)
#   Pos 133-134 : Population year (2 chars)
#   Pos 135-141 : Enrollment (7 chars, numeric)
#   Pos 142-143 : Enrollment year (2 chars)
#   Pos 144-145 : Function code for special districts (2 chars)
#   Pos 146-147 : School level code (2 chars)
#   Pos 148-151 : Fiscal year ending (4 chars)
#   Pos 152-153 : Survey year (2 chars)
# =============================================================================

library(readr)
library(dplyr)
library(stringr)
library(here)

# -----------------------------------------------------------------------------
# 1.  File path – update this to wherever Fin_GID_2012.txt lives on your system
# -----------------------------------------------------------------------------

fin_path <- here("DATA", "2012 COG-F", "2012FinEstDAT_07242018modp_pu.txt")

# --- State code lookup -------------------------------------------------------
# NOTE: These are COG-F *internal* state codes — they are NOT Census FIPS codes.
# The COG-F numbering differs from Census FIPS for many states, e.g.:
#   Idaho:      COG-F = "13",  Census FIPS = "16"
#   Montana:    COG-F = "27",  Census FIPS = "30"
#   Wyoming:    COG-F = "51",  Census FIPS = "56"
# Use these codes only for filtering the finance file (state_code field).
# Census FIPS state codes come from the GID file (fips_state, positions 114-115)
# and are the correct key for joining to external datasets (PAD-US, ACS, TIGER).
state_codes <- c(
  "00"="United States", "01"="Alabama",        "02"="Alaska",
  "03"="Arizona",       "04"="Arkansas",       "05"="California",
  "06"="Colorado",      "07"="Connecticut",    "08"="Delaware",
  "09"="District of Columbia", "10"="Florida", "11"="Georgia",
  "12"="Hawaii",        "13"="Idaho",          "14"="Illinois",
  "15"="Indiana",       "16"="Iowa",           "17"="Kansas",
  "18"="Kentucky",      "19"="Louisiana",      "20"="Maine",
  "21"="Maryland",      "22"="Massachusetts",  "23"="Michigan",
  "24"="Minnesota",     "25"="Mississippi",    "26"="Missouri",
  "27"="Montana",       "28"="Nebraska",       "29"="Nevada",
  "30"="New Hampshire", "31"="New Jersey",     "32"="New Mexico",
  "33"="New York",      "34"="North Carolina", "35"="North Dakota",
  "36"="Ohio",          "37"="Oklahoma",       "38"="Oregon",
  "39"="Pennsylvania",  "40"="Rhode Island",   "41"="South Carolina",
  "42"="South Dakota",  "43"="Tennessee",      "44"="Texas",
  "45"="Utah",          "46"="Vermont",        "47"="Virginia",
  "48"="Washington",    "49"="West Virginia",  "50"="Wisconsin",
  "51"="Wyoming"
)

# --- Type-of-government lookup -----------------------------------------------
gov_types <- c(
  "0" = "State",
  "1" = "County",
  "2" = "City",
  "3" = "Township",
  "4" = "Special District",
  "5" = "Independent School District or Educational Service Agency"
)

# --- Item code descriptions (selected key codes) -----------------------------
# Full list in technical documentation; a representative sample is included here.
item_descriptions <- c(
  "T01"="Property Tax",
  "T09"="Total General Sales Taxes",
  "T40"="Indiv Income Tax",
  "T41"="Corporation Net Income Tax",
  "E12"="Elem Sec Educ-Current Oper",
  "E18"="Higher Educ Current Operation",
  "E36"="Hospitals-Current Oper",
  "E44"="Regular Highways-Current Oper",
  "E62"="Police Protection-Current Oper",
  "E80"="Sewerage-Current Operation",
  "A91"="Water Utilities Revenue",
  "A92"="Electric Utilities Revenue",
  "B21"="Federal IG Education",
  "B46"="Federal IG Highways",
  "B79"="Federal IG Public Welfare",
  "C21"="State IG Education",
  "Z00"="Total Salaries & Wages",
  "I89"="General Interest on Debt"
  # Add more as needed from the technical documentation
)

# =============================================================================
# PARSE FUNCTION
# =============================================================================

parse_fin_est <- function(fin_path,
                          state_filter   = NULL,   # western states
                          govtype_filter = NULL,   # all except state
                          item_filter    = NULL,   # e.g. c("T01","T09")
                          add_labels     = TRUE) {

  message("Reading file: ", fin_path)
  message("File size: ", round(file.size(fin_path) / 1e6, 1), " MB")

  # Read all lines as raw fixed-width text
  raw <- read_lines(fin_path, progress = TRUE)
  message("Total records read: ", format(length(raw), big.mark = ","))

  # Validate record length
  rec_lengths <- nchar(raw)
  bad <- sum(rec_lengths != 34)
  if (bad > 0) warning(bad, " records do not have expected length of 34 characters.")

  # Extract fields using substr (positions are 1-based per documentation)
  # IMPORTANT: state_code (positions 1-2) is the COG-F *internal* state code,
  # NOT a Census FIPS code. Idaho is "13" here but "16" in Census FIPS.
  # Census FIPS codes are stored separately in the GID file (fips_state,
  # positions 114-115) and are used for all external joins (PAD-US, ACS, TIGER).
  df <- tibble(
    id            = substr(raw,  1, 14),
    state_code    = substr(raw,  1,  2),  # COG-F internal code — NOT Census FIPS
    gov_type_code = substr(raw,  3,  3),
    county_code   = substr(raw,  4,  6),
    unit_id       = substr(raw,  7,  9),
    unit_suffix   = substr(raw, 10, 14),   # should be "00000"
    item_code     = substr(raw, 15, 17),
    amount_raw    = substr(raw, 18, 29),
    year          = substr(raw, 30, 33),
    flag          = substr(raw, 34, 34)
  )

  # Convert amount to numeric (stored in thousands of dollars)
  df <- df |>
    mutate(amount_thousands = as.numeric(trimws(amount_raw))) |>
    select(-amount_raw)

  # Apply filters before adding labels (faster on large file)
  if (!is.null(state_filter)) {
    df <- filter(df, state_code %in% state_filter)
    message("After state filter: ", format(nrow(df), big.mark = ","), " records")
  }
  if (!is.null(govtype_filter)) {
    df <- filter(df, gov_type_code %in% govtype_filter)
    message("After gov-type filter: ", format(nrow(df), big.mark = ","), " records")
  }
  if (!is.null(item_filter)) {
    df <- filter(df, item_code %in% item_filter)
    message("After item filter: ", format(nrow(df), big.mark = ","), " records")
  }

  # Optionally add human-readable labels
  if(add_labels) {
    df <- df |>
      mutate(
        state_name    = state_codes[state_code],
        gov_type_name = gov_types[gov_type_code],
        item_desc     = item_descriptions[item_code],   # NA if not in lookup
        flag_desc     = case_when(
          flag == "I" ~ "Imputed",
          flag == "R" ~ "Reported",
          TRUE        ~ paste0("Unknown (", flag, ")"))
      )
  }

  message("Final record count: ", format(nrow(df), big.mark = ","))
  df
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

fin_path <- "2012FinEstDAT_07242018modp_pu.txt"

# --- Example 1: Load everything (all states, all gov types, all items) -------
# df_all <- parse_fin_est(fin_path)

# --- Example 2: State-level records only (gov_type = "0") --------------------
df_states <- parse_fin_est(
  fin_path,
  govtype_filter = "0"   # State governments only
)

# --- Example 3: Property and sales tax for Colorado (state code "06") --------
df_co_tax <- parse_fin_est(
  fin_path,
  state_filter = "06",
  item_filter  = c("T01", "T09", "T40")
)

# --- Example 4: Education spending across all governments --------------------
df_education <- parse_fin_est(
  fin_path,
  item_filter = c("E12", "E18", "B21", "C21")
)

# =============================================================================
# BASIC SUMMARIES
# =============================================================================

# Total property tax revenue by state (amounts in $thousands)
property_tax_by_state <- df_states |>
  filter(item_code == "T01") |>
  arrange(desc(amount_thousands)) |>
  select(state_name, state_code, item_code, amount_thousands)

print("Top 10 states by property tax revenue ($thousands):")
print(head(property_tax_by_state, 10))

# Summary of Colorado tax data
if (nrow(df_co_tax) > 0) {
  co_summary <- df_co_tax |>
    group_by(item_code, item_desc) |>
    summarise(
      total_thousands = sum(amount_thousands, na.rm = TRUE),
      n_govts         = n_distinct(id),
      .groups = "drop"
    )
  print("Colorado tax summary:")
  print(co_summary)
}

# =============================================================================
# EXPORT (optional)
# =============================================================================

# Write a filtered subset to CSV
# write_csv(df_co_tax, "2012_colorado_tax.csv")
# write_csv(df_states, "2012_state_governments.csv")

# =============================================================================
# GID PARSE FUNCTION
# =============================================================================

parse_gid <- function(gid_path) {

  message("Reading GID file: ", gid_path)
  message("File size: ", round(file.size(gid_path) / 1e6, 1), " MB")

  # Read lines; strip Windows \r\n carriage returns
  raw <- read_lines(gid_path, progress = TRUE)
  raw <- gsub("\r", "", raw, fixed = TRUE)

  message("Total GID records read: ", format(length(raw), big.mark = ","))

  # Validate record length (expected: 153 chars after stripping \r)
  rec_lengths <- nchar(raw)
  bad <- sum(rec_lengths != 153)
  if (bad > 0) warning(bad, " GID records do not have expected length of 153 characters.")

  gid <- tibble(
    id                 = substr(raw,   1,  14),
    govt_name          = trimws(substr(raw,  15,  78)),
    county_name        = trimws(substr(raw,  79, 113)),
    # fips_state and fips_county are Census FIPS codes (positions 114-118).
    # These differ from the COG-F internal state_code in the finance file —
    # e.g. Idaho's Census FIPS is "16" but its COG-F state_code is "13".
    # fips_state is the correct key for joining to PAD-US, ACS, and TIGER data.
    fips_state         = substr(raw, 114, 115),  # Census FIPS state code
    fips_county        = substr(raw, 116, 118),  # Census FIPS county code
    fips_place         = substr(raw, 119, 123),
    population         = as.numeric(trimws(substr(raw, 124, 132))),
    population_year    = trimws(substr(raw, 133, 134)),
    enrollment         = as.numeric(trimws(substr(raw, 135, 141))),
    enrollment_year    = trimws(substr(raw, 142, 143)),
    special_dist_func  = trimws(substr(raw, 144, 145)),
    school_level_code  = trimws(substr(raw, 146, 147)),
    fiscal_year_ending = trimws(substr(raw, 148, 151)),
    survey_year        = trimws(substr(raw, 152, 153))
  )

  message("GID parse complete: ", format(nrow(gid), big.mark = ","), " records")
  gid
}

# =============================================================================
# JOINED DATASET FUNCTION
# =============================================================================
# Left-joins the finance file (one row per govt x item code) to the GID
# directory (one row per govt) on the 14-character ID key, enriching every
# finance record with govt name, county, population, enrollment, and fiscal
# year metadata.

parse_and_join <- function(fin_path,
                           gid_path,
                           state_filter   = NULL,
                           govtype_filter = NULL,
                           item_filter    = NULL,
                           add_labels     = TRUE) {

  fin <- parse_fin_est(fin_path,
                       state_filter   = state_filter,
                       govtype_filter = govtype_filter,
                       item_filter    = item_filter,
                       add_labels     = add_labels)

  gid <- parse_gid(gid_path)

  # Left join: retain all finance rows; attach GID columns where matched
  joined <- fin |>
    left_join(
      gid |> select(id, govt_name, county_name, fips_state, fips_county,
                    fips_place, population, population_year, enrollment,
                    enrollment_year, special_dist_func, school_level_code,
                    fiscal_year_ending, survey_year),
      by = "id"
    )

  unmatched <- sum(is.na(joined$govt_name))
  if (unmatched > 0) {
    message("Note: ", format(unmatched, big.mark = ","),
            " finance records had no matching GID entry.")
  }

  message("Joined dataset: ", format(nrow(joined), big.mark = ","),
          " rows x ", ncol(joined), " columns")
  joined
}

# =============================================================================
# JOINED EXAMPLES
# =============================================================================

fin_path <- "2012FinEstDAT_07242018modp_pu.txt"
gid_path <- "Fin_GID_2012.txt"

# --- Example A: All state governments, all items, with GID metadata ----------
df_states_joined <- parse_and_join(
  fin_path, gid_path,
  govtype_filter = "0"
)

# --- Example B: Idaho only, all gov types, key tax items ------------------
df_id_joined <- parse_and_join(
  fin_path, gid_path,
  state_filter = "13",
  item_filter  = c("T01", "T09", "T40", "E12", "E62")
)

# --- Example C: Police spending for cities, with population for per-capita ---
df_police_cities <- parse_and_join(
  fin_path, gid_path,
  govtype_filter = "2",          # Cities only
  item_filter    = "E62"         # Police Protection - Current Oper
) |>
  mutate(
    spend_per_capita = ifelse(
      !is.na(population) & population > 0,
      (amount_thousands * 1000) / population,
      NA_real_
    )
  )

# =============================================================================
# JOINED SUMMARIES
# =============================================================================

# Top 10 cities by police spending per capita
print("Top 10 cities by police spending per capita ($):")
df_police_cities |>
  filter(!is.na(spend_per_capita)) |>
  arrange(desc(spend_per_capita)) |>
  select(govt_name, state_name, county_name, population,
         amount_thousands, spend_per_capita) |>
  head(10) |>
  print()

# Idaho joined summary: total spending by item and govt type
print("Idaho: total by item code and govt type:")
df_id_joined |>
  group_by(item_code, item_desc, gov_type_name) |>
  summarise(
    total_thousands = sum(amount_thousands, na.rm = TRUE),
    n_govts         = n_distinct(id),
    .groups = "drop"
  ) |>
  arrange(item_code, gov_type_name) |>
  print()

# =============================================================================
# EXPORT (optional)
# =============================================================================

# write_csv(df_states_joined,  "2012_states_joined.csv")
# write_csv(df_co_joined,      "2012_colorado_joined.csv")
# write_csv(df_police_cities,  "2012_police_cities_joined.csv")

# =============================================================================
# WESTERN STATES FILTER: cogf_2012w
# States: AZ, CA, CO, ID, MT, NM, NV, OR, UT, WA, WY
# =============================================================================

western_state_codes <- c(
  "03",  # Arizona       (AZ)
  "05",  # California    (CA)
  "06",  # Colorado      (CO)
  "13",  # Idaho         (ID)
  "27",  # Montana       (MT)
  "32",  # New Mexico    (NM)
  "29",  # Nevada        (NV)
  "38",  # Oregon        (OR)
  "45",  # Utah          (UT)
  "48",  # Washington    (WA)
  "51"   # Wyoming       (WY)
)

message("Building cogf_2012w: western states combined dataset...")

cogf_2012w <- parse_and_join(
  fin_path     = here("DATA", "2012 COG-F", "2012FinEstDAT_07242018modp_pu.txt"),
  gid_path     = here("DATA", "2012 COG-F", "Fin_GID_2012.txt"),
  state_filter = western_state_codes,
  add_labels   = TRUE
) |>
  filter(gov_type_code != "0")  # Remove state-level governments

message("cogf_2012w: ", format(nrow(cogf_2012w), big.mark = ","),
        " rows x ", ncol(cogf_2012w), " columns")

# Quick summary by state
cogf_2012w |>
  group_by(state_code, state_name) |>
  summarise(
    n_records   = n(),
    n_govts     = n_distinct(id),
    n_items     = n_distinct(item_code),
    .groups = "drop"
  ) |>
  arrange(state_name) |>
  print()

# Add class_type based on item_code prefix (lookup: Class Types.txt)
cogf_2012w <- cogf_2012w |>
  mutate(
    class_type = case_when(
      substr(item_code, 1, 1) == "T"                                  ~ "R_Tax",
      substr(item_code, 1, 1) %in% c("B", "C", "D")                  ~ "R_IG",
      substr(item_code, 1, 1) == "A"                                  ~ "R_Ch",
      substr(item_code, 1, 1) == "U"                                  ~ "R_Oth",
      substr(item_code, 1, 1) %in% c("E", "I", "J", "F", "G", "K")  ~ "E_Dir",
      substr(item_code, 1, 1) %in% c("L", "M", "Q", "S")             ~ "E_IG",
      grepl("^[0-9]", item_code) ~ "D"  # Codes beginning with a number: Debt
    )
  )

# Save as .rds (preserves column types; fast to reload with readRDS())
saveRDS(cogf_2012w, here("R Code and Data", "cogf_2012w.rds"))
message("Saved: cogf_2012w.rds")

# Also save as CSV if a plain-text copy is needed
# write_csv(cogf_2012w, here("R Code and Data", "2012_All_Gov_West.csv"))

# Save as xlsx
library(openxlsx)
write.xlsx(cogf_2012w, file = here("R Code and Data", "2012_All_Gov_West.xlsx"), sheetName = "2012_All_Gov_West", rowNames = FALSE)
message("Saved: 2012_All_Gov_West.xlsx")
