# 07_plots.R ==================================================================
# Purpose: Correlations between federal land share and per-capita COG-F
#          spending by item code and government type.
#          Produces bar charts of top-10 correlated item codes for:
#            (1) Full western-states sample, by government type
#            (2) Non-metro counties only, by government type
#          Also produces scatter plots of land share vs. county-level tax
#          and expenditure fields from Cty_Agg.rds.
#
# Inputs:  cogf_west.rds (2012 wave), fed_share_west.csv, ruralurbancodes2013.xls,
#          Cty_Agg.rds
# =============================================================================

library(tidyverse)
library(here)
library(readxl)
library(tidytext)   # reorder_within / scale_y_reordered

# --- Item code descriptions (key codes; extend as needed) --------------------
item_descriptions <- c(
  "T01" = "Property Tax",
  "T09" = "Total General Sales Taxes",
  "T40" = "Indiv Income Tax",
  "T41" = "Corp Net Income Tax",
  "E12" = "Elem Sec Educ-Current Oper",
  "E18" = "Higher Educ Current Oper",
  "E36" = "Hospitals-Current Oper",
  "E44" = "Regular Highways-Current Oper",
  "E62" = "Police Protection-Current Oper",
  "E80" = "Sewerage-Current Oper",
  "A91" = "Water Utilities Revenue",
  "A92" = "Electric Utilities Revenue",
  "B21" = "Federal IG Education",
  "B46" = "Federal IG Highways",
  "B79" = "Federal IG Public Welfare",
  "C21" = "State IG Education",
  "Z00" = "Total Salaries & Wages",
  "I89" = "General Interest on Debt"
)

# --- State FIPS lookup (western states only) ---------------------------------
state_fips_lookup <- tibble(
  state_name = c("Arizona",   "California", "Colorado", "Idaho",    "Montana",
                 "Nevada",    "New Mexico", "Oregon",   "Utah",     "Washington", "Wyoming"),
  state_fips = c("04",        "06",         "08",       "16",       "30",
                 "32",        "35",         "41",       "49",       "53",         "56")
)

# =============================================================================
# 1.  LOAD & JOIN DATA
# =============================================================================

cogf <- readRDS(here("output", "cogf_west.rds")) |> filter(year == 2012)

# Federal land share (%)
fed_share <- read.csv(
  here("output", "fed_share_west.csv"),
  colClasses = c(GEOID = "character")
) |>
  select(GEOID, fed_share_pct, usfs_share_pct, usfs_doi_share_pct)

# ERS Rural-Urban Continuum Codes 2013
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
  select(GEOID, metro)

# Join fed share and metro status; shorten ISD label for plots
cogf_fed <- cogf |>
  left_join(state_fips_lookup, by = "state_name") |>
  mutate(GEOID = if_else(!is.na(state_fips), paste0(state_fips, fips_county), NA_character_)) |>
  left_join(fed_share,   by = "GEOID") |>
  left_join(metro_join,  by = "GEOID") |>
  select(-state_fips, -GEOID) |>
  mutate(gov_type = case_when(
    gov_type_name == "Independent School District or Educational Service Agency" ~
      "Independent\nSchool District",
    TRUE ~ gov_type_name
  ))

# County-level population (from county government records)
county_pop <- cogf_fed |>
  filter(gov_type_code == "1") |>
  distinct(fips_state, fips_county, population)

# County-level federal land share
county_fed_df <- cogf_fed |>
  distinct(fips_state, fips_county, fed_share_pct, usfs_share_pct, usfs_doi_share_pct)

# Non-metro subset and its county-level reference tables
cogf_nonmetro <- cogf_fed |> filter(metro == "Non-Metro")

county_pop_nm <- county_pop |>
  semi_join(
    cogf_nonmetro |> distinct(fips_state, fips_county),
    by = c("fips_state", "fips_county")
  )

county_fed_nm <- county_fed_df |>
  semi_join(
    cogf_nonmetro |> distinct(fips_state, fips_county),
    by = c("fips_state", "fips_county")
  )

# =============================================================================
# 2.  COMPUTE CORRELATIONS
#     For each gov_type x item_code: aggregate to county level, compute
#     per-capita spending, then correlate with federal land share (%).
# =============================================================================

share_labels <- c(
  fed_share_pct      = "Federal Land Share",
  usfs_share_pct     = "USFS Land Share",
  usfs_doi_share_pct = "USFS + DOI Land Share"
)

compute_cors <- function(data, pop_df, fed_df) {
  data |>
    group_by(fips_state, fips_county, gov_type, item_code) |>
    summarise(amount = sum(amount_thousands, na.rm = TRUE), .groups = "drop") |>
    left_join(pop_df, by = c("fips_state", "fips_county")) |>
    left_join(fed_df, by = c("fips_state", "fips_county")) |>
    mutate(amount_pc = amount / population) |>
    filter(!is.na(amount_pc), is.finite(amount_pc), population > 0) |>
    pivot_longer(
      cols      = c(fed_share_pct, usfs_share_pct, usfs_doi_share_pct),
      names_to  = "share_type",
      values_to = "share_pct"
    ) |>
    filter(!is.na(share_pct)) |>
    group_by(share_type, gov_type, item_code) |>
    summarise(
      # NA (not a warning) when either variable is constant across counties
      r = if (isTRUE(sd(share_pct) > 0 && sd(amount_pc) > 0)) cor(share_pct, amount_pc) else NA_real_,
      n = n(),
      .groups = "drop"
    ) |>
    mutate(item_desc = item_descriptions[item_code])
}

cors_by_type  <- compute_cors(cogf_fed,      county_pop,    county_fed_df)
cors_nonmetro <- compute_cors(cogf_nonmetro, county_pop_nm, county_fed_nm)

# =============================================================================
# 3.  SELECT TOP 10 PER GOVERNMENT TYPE
#     Exclude aggregate/rollup item codes (numeric suffix 90-99).
# =============================================================================

build_top10 <- function(cors_df) {
  cors_df |>
    filter(!grepl("9[0-9]$", item_code), !is.na(r)) |>
    group_by(share_type, gov_type) |>
    slice_max(abs(r), n = 10) |>
    arrange(share_type, gov_type, desc(abs(r))) |>
    mutate(label = if_else(
      !is.na(item_desc),
      paste0(item_code, ": ", str_trunc(item_desc, 20)),
      item_code
    )) |>
    ungroup()
}

top_by_type  <- build_top10(cors_by_type)
top_nonmetro <- build_top10(cors_nonmetro)

# =============================================================================
# 4.  BAR CHARTS
# =============================================================================

plot_top_cors <- function(top_df, stype, title) {
  x_label <- paste0("Pearson r  (", share_labels[stype], " % vs. per-capita spending)")
  top_df |>
    filter(share_type == stype) |>
    mutate(label = reorder_within(label, r, gov_type)) |>
    ggplot(aes(x = r, y = label)) +
    geom_col() +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    scale_y_reordered() +
    facet_wrap(~gov_type, scales = "free_y") +
    labs(title = title, x = x_label, y = NULL)
}

# All western counties
plot_top_cors(top_by_type,  "fed_share_pct",
  "Top 10 COG-F Item Codes by Government Type\nAll Western Counties — Federal Land Share")

plot_top_cors(top_by_type,  "usfs_share_pct",
  "Top 10 COG-F Item Codes by Government Type\nAll Western Counties — USFS Land Share")

plot_top_cors(top_by_type,  "usfs_doi_share_pct",
  "Top 10 COG-F Item Codes by Government Type\nAll Western Counties — USFS + DOI Land Share")

# Non-metro counties only
plot_top_cors(top_nonmetro, "fed_share_pct",
  "Top 10 COG-F Item Codes by Government Type\nNon-Metro Counties Only — Federal Land Share")

plot_top_cors(top_nonmetro, "usfs_share_pct",
  "Top 10 COG-F Item Codes by Government Type\nNon-Metro Counties Only — USFS Land Share")

plot_top_cors(top_nonmetro, "usfs_doi_share_pct",
  "Top 10 COG-F Item Codes by Government Type\nNon-Metro Counties Only — USFS + DOI Land Share")

# =============================================================================
# 5.  SCATTER PLOTS: LAND SHARE vs. TAX & EXPENDITURE FIELDS (Cty_Agg)
# =============================================================================

cty_agg <- readRDS(here("output", "Cty_Agg.rds")) |>
  filter(year == 2012) |>
  # Per-capita versions of raw dollar fields ($000s per person)
  mutate(across(
    c(R_Tax, R_IG, R_IG_Fed, R_Ch, R_Oth, E_Dir, E_IG, Total_Revenue, Total_Expenditures),
    ~ . / population,
    .names = "PerCap_{.col}"
  ))

fiscal_labels <- c(
  PerCap_R_Tax              = "PC Tax Revenue",
  PerCap_R_IG               = "PC IGR Revenue",
  PerCap_R_IG_Fed           = "PC Federal IG Revenue",
  PerCap_R_Ch               = "PC Charges Revenue",
  PerCap_R_Oth              = "PC Misc Revenue",
  PerCap_E_Dir              = "PC Direct Expenditures",
  PerCap_E_IG               = "PC IG Expenditures",
  PerCap_Total_Revenue      = "PC Total Revenue",
  PerCap_Total_Expenditures = "PC Total Expenditures",
  PerCap_TaxBurden          = "PC Tax Burden (Tax+Charges)",
  PerCap_Expenditures       = "PC Expenditures (Direct)"
)

# Pivot share types and fiscal variables to long form
# any_of() on share cols allows graceful handling if usfs_doi_share_pct is not
# yet present (requires re-running 02_fed_share_acs.R and 04_Cty_Agg.R first)
cty_long <- cty_agg |>
  pivot_longer(
    cols      = any_of(names(share_labels)),
    names_to  = "share_type",
    values_to = "share_pct"
  ) |>
  pivot_longer(
    cols      = all_of(names(fiscal_labels)),
    names_to  = "fiscal_var",
    values_to = "fiscal_val"
  ) |>
  mutate(
    share_label  = factor(share_labels[share_type], levels = share_labels),
    fiscal_label = factor(fiscal_labels[fiscal_var], levels = fiscal_labels)
  ) |>
  filter(!is.na(share_pct), !is.na(fiscal_val))

# One plot per fiscal variable; panels for each share type
walk(levels(cty_long$fiscal_label), function(fvar) {
  p <- cty_long |>
    filter(fiscal_label == fvar) |>
    ggplot(aes(x = share_pct, y = fiscal_val)) +
    geom_point(alpha = 0.35, size = 0.9) +
    facet_wrap(~share_label, scales = "free_x", ncol = 3) +
    labs(
      title = fvar,
      x     = "Land Share (%)",
      y     = fvar
    )
  print(p)
})

# =============================================================================
# 6.  NON-METRO DEEP DIVE: PC Tax Revenue vs. USFS + DOI Land Share
# =============================================================================

# Cty_Agg already carries usfs_doi_share_pct and metro
cty_agg |>
  dplyr::filter(metro == "Non-Metro", !is.na(usfs_doi_share_pct), !is.na(PerCap_R_Tax)) |>
  ggplot(aes(x = usfs_doi_share_pct, y = PerCap_R_Tax)) +
  geom_point(alpha = 0.4) +
  labs(
    title = "Per Capita Tax Revenue vs. USFS + DOI Land Share\nNon-Metro Counties Only",
    x     = "USFS + DOI Land Share (%)",
    y     = "PC Tax Revenue ($000s per person)"
  )
