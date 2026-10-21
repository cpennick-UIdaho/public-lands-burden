
library(tidyverse)
library(openxlsx)

data1 <- readRDS("cogf_2012w.rds")

# 1. Count of distinct governments per county
gov_counts <- data1 |>
  distinct(state_code, county_name, id) |>
  count(state_code, county_name, name = "n_govts")

# 2. Aggregate amount_thousands by class_type (Class Code) and county
cty_class_agg <- data1 |>
  filter(!is.na(class_type)) |>
  group_by(state_code, county_name, class_type) |>
  summarise(total_thousands = sum(amount_thousands, na.rm = TRUE), .groups = "drop")

# 3. Join government counts and state_name onto aggregated data
state_names <- data1 |>
  distinct(state_code, state_name)

data2 <- cty_class_agg |>
  left_join(gov_counts, by = c("state_code", "county_name")) |>
  left_join(state_names, by = "state_code") |>
  pivot_wider(names_from = class_type, values_from = total_thousands) |>
  select(state_code, state_name, county_name, n_govts, R_Tax, R_IG, R_Ch, R_X, E_Dir, E_IG) |>
  mutate(
    Total_Revenue      = rowSums(across(c(R_Tax, R_IG, R_Ch, R_X)), na.rm = TRUE),
    Total_Expenditures = rowSums(across(c(E_Dir, E_IG)), na.rm = TRUE)
  )

# Save as .rds
saveRDS(data2, "Cty_Agg.rds")
message("Saved: Cty_Agg.rds")

# Also save as CSV if a plain-text copy is needed
# write_csv(data2, "County_Aggregate.csv")

# Save as xlsx with full column names for readability
data2_xlsx <- data2 |>
  rename(
    "General Revenue-Taxes"   = R_Tax,
    "General Revenue-IGR"     = R_IG,
    "General Revenue-Charges" = R_Ch,
    "General Revenue-Misc"    = R_X,
    "Expenditures-Direct"     = E_Dir,
    "Expenditures-IG"         = E_IG,
    "Total Revenue"           = Total_Revenue,
    "Total Expenditures"      = Total_Expenditures
  )

write.xlsx(data2_xlsx, file = "County_Aggregate.xlsx", sheetName = "Cty_Agg", rowNames = FALSE)
message("Saved: County_Aggregate.xlsx")
