# =============================================================================
# Parse Fin_GID_2012.txt
# 2012 State & Local Government Finances – GID Directory Information File
# Source: U.S. Census Bureau
#
# File format: ASCII fixed-width, 153 characters per record (no delimiter)
# Technical documentation: 2012_S_L_Indiv_Unit_Data_File_Tech_Doc.pdf
# =============================================================================

library(readr)     # read_fwf()
library(dplyr)     # data wrangling
library(stringr)   # string trimming helpers
library(openxlsx)  # write formatted .xlsx output

# -----------------------------------------------------------------------------
# 1.  File path – update this to wherever Fin_GID_2012.txt lives on your system
# -----------------------------------------------------------------------------

setwd("C:/Users/chelseap/OneDrive - University of Idaho/Projects/NSF ERF/COG-F Data/2012 COG-F")
gid_2012 <- "Fin_GID_2012.txt"

# -----------------------------------------------------------------------------
# 2.  Fixed-width column specification (positions are 1-based, per the tech doc)
#
#  Field                              Start  End  Width  Type
#  ---------------------------------- -----  ---  -----  ----
#  1  ID code                             1   14    14   Char
#  2  ID name                            15   78    64   Char
#  3  County name                        79  113    35   Char
#  4  FIPS state code                   114  115     2   Char
#  5  FIPS county code                  116  118     3   Char
#  6  FIPS place code                   119  123     5   Char
#  7  Population                        124  132     9   Num
#  8  Population year                   133  134     2   Char
#  9  Enrollment                        135  141     7   Num
# 10  Enrollment year                   142  143     2   Char
# 11  Function code (special districts) 144  145     2   Char
# 12  School level code                 146  147     2   Char
# 13  Fiscal year ending                148  151     4   Char
# 14  Survey year                       152  153     2   Char
# -----------------------------------------------------------------------------

col_positions <- fwf_cols(
  fips_state       = c(1, 2),
  county_name      = c(79,  113),
  gov_type         = c(3, 3),
  fips_county      = c(4, 6),
  unit_id          = c(7, 12),
  fips_place       = c(119, 123),
  population       = c(124, 132),
  population_year  = c(133, 134),
  enrollment       = c(135, 141),
  enrollment_year  = c(142, 143),
  function_code    = c(144, 145),
  school_level     = c(146, 147),
  fiscal_year_end  = c(148, 151),
  survey_year      = c(152, 153)
)

col_types <- cols(
  id_code         = col_character(),
  id_name         = col_character(),
  county_name     = col_character(),
  fips_state      = col_character(),
  fips_county     = col_character(),
  fips_place      = col_character(),
  population      = col_double(),    # 9-digit field; double avoids integer overflow
  population_year = col_character(),
  enrollment      = col_double(),
  enrollment_year = col_character(),
  function_code   = col_character(),
  school_level    = col_character(),
  fiscal_year_end = col_character(),
  survey_year     = col_character()
)

# -----------------------------------------------------------------------------
# 3.  Read the file
# -----------------------------------------------------------------------------
message("Reading ", gid_2012, " ...")

gid <- read_fwf(
  file      = gid_2012,
  col_positions = col_positions,
  col_types = col_types,
  trim_ws   = TRUE   # strip leading/trailing whitespace from all fields
)

message("Rows read: ", nrow(gid))

# -----------------------------------------------------------------------------
# 4.  Parse the composite ID code into its component parts
#
#  Positions within id_code (1-based within the 14-char field):
#    1-2   State code (numeric, matches FIPS state code)
#    3     Type-of-government code
#            0 = State
#            1 = County
#            2 = City
#            3 = Township
#            4 = Special District
#            5 = Independent School District / Educational Service Agency
#    4-6   County or county-type area code
#    7-9   Unit identifier
#   10-14  Should be "00000" (unit is not part of another government)
# -----------------------------------------------------------------------------

gov_type_labels <- c(
  "0" = "State",
  "1" = "County",
  "2" = "City",
  "3" = "Township",
  "4" = "Special District",
  "5" = "Independent School District / ESA"
)

gid <- gid |>
  mutate(
    id_state_code  = substr(id_code, 1,  2),
    id_gov_type    = substr(id_code, 3,  3),
    id_county_code = substr(id_code, 4,  6),
    id_unit_id     = substr(id_code, 7,  9),
    id_suffix      = substr(id_code, 10, 14),
    gov_type_label = gov_type_labels[id_gov_type]
  )

# -----------------------------------------------------------------------------
# 5.  Recode blank / zero numeric fields to NA
#     Blank spaces read as NA by read_fwf; zero population/enrollment may be
#     legitimate, but the file uses spaces for "not applicable" fields.
# -----------------------------------------------------------------------------

gid <- gid |>
  mutate(
    population  = if_else(population  == 0 & is.na(population_year),  NA_real_, population),
    enrollment  = if_else(enrollment  == 0 & is.na(enrollment_year),  NA_real_, enrollment)
  )

# -----------------------------------------------------------------------------
# 6.  Quick validation
# -----------------------------------------------------------------------------

cat("\n--- Structure ---\n")
glimpse(gid)

cat("\n--- Government type distribution ---\n")
print(count(gid, gov_type_label, sort = TRUE))

cat("\n--- FIPS state code distribution (top 10) ---\n")
print(count(gid, fips_state, sort = TRUE) |> head(10))

cat("\n--- First 5 rows ---\n")
print(head(gid, 5))

# -----------------------------------------------------------------------------
# 7.  Save parsed data as a formatted Excel workbook
# -----------------------------------------------------------------------------

output_xlsx <- sub("\\.txt$", "_parsed.xlsx", gid_2012)

wb <- createWorkbook()
addWorksheet(wb, "COG-F 2012")

# --- Styles ------------------------------------------------------------------
header_style <- createStyle(
  fontName      = "Arial",
  fontSize      = 10,
  fontColour    = "white",
  fgFill        = "#2F5496",        # dark blue header background
  halign        = "center",
  valign        = "center",
  textDecoration = "bold",
  border        = "Bottom",
  borderColour  = "#FFFFFF",
  wrapText      = TRUE
)

body_style <- createStyle(
  fontName = "Arial",
  fontSize = 10
)

number_style <- createStyle(
  fontName   = "Arial",
  fontSize   = 10,
  numFmt     = "#,##0"              # comma-separated thousands
)

# --- Write data --------------------------------------------------------------
writeData(wb, "COG-F 2012", gid, startRow = 1, startCol = 1, headerStyle = header_style)

# Apply body style to all data rows
addStyle(wb, "COG-F 2012", body_style,
         rows    = 2:(nrow(gid) + 1),
         cols    = 1:ncol(gid),
         gridExpand = TRUE)

# Apply number style to numeric columns (population, enrollment)
num_cols <- which(names(gid) %in% c("population", "enrollment"))
addStyle(wb, "COG-F 2012", number_style,
         rows    = 2:(nrow(gid) + 1),
         cols    = num_cols,
         gridExpand = TRUE)

# --- Column widths -----------------------------------------------------------
col_widths <- c(
  id_code         = 16,
  id_name         = 40,
  county_name     = 25,
  fips_state      = 12,
  fips_county     = 14,
  fips_place      = 14,
  population      = 14,
  population_year = 16,
  enrollment      = 14,
  enrollment_year = 16,
  function_code   = 16,
  school_level    = 14,
  fiscal_year_end = 16,
  survey_year     = 12,
  id_state_code   = 14,
  id_gov_type     = 14,
  id_county_code  = 16,
  id_unit_id      = 12,
  id_suffix       = 12,
  gov_type_label  = 36
)

for (col_name in names(col_widths)) {
  col_idx <- which(names(gid) == col_name)
  if (length(col_idx) > 0) {
    setColWidths(wb, "COG-F 2012", cols = col_idx, widths = col_widths[col_name])
  }
}

# Freeze the top header row
freezePane(wb, "COG-F 2012", firstRow = TRUE)

# --- Save --------------------------------------------------------------------
saveWorkbook(wb, output_xlsx, overwrite = TRUE)
message("\nParsed data written to: ", output_xlsx)
