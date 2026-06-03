# ======== USER INPUT ========
input_xlsx <- "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/1. info/Demographics_EPL_unprocessed_removed.xlsx"
output_png <- "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/1. info/table1.png"
# ============================

# --- Inspect Case column to diagnose ---
cat("\nUnique Case values:\n")
print(unique(df$Case))

cat("\nRows where Case is missing or non-numeric:\n")
print(df %>% filter(is.na(Case) | !(Case %in% c(0,1))) %>% select(StudyID, Case))

nrow(df)

# --- Libraries ---
library(readxl)
library(dplyr)
library(gt)
library(webshot2)   # install.packages("webshot2") if not installed

# --- Read & clean ---
# df <- read_excel(input_xlsx)
df <- read_excel(input_xlsx, sheet = 2)

names(df) <- trimws(names(df))

df <- df %>%
  mutate(
    across(c(Case, Parous, ConceptionMode, Race_Ethnicity), ~as.numeric(.)),
    Case = ifelse(is.na(Case), 0, Case)
  )

# --- Split groups ---
cases <- df %>% filter(Case == 1)
controls <- df %>% filter(Case == 0)

cat("\nCases:", nrow(cases), "Controls:", nrow(controls), "\n")

# --- Helper to summarize group ---
summarize_group <- function(dat, group_name) {
  n <- nrow(dat)
  out <- list(
    "Age (year), Mean ± SD" = sprintf("%.2f ± %.2f",
                                      mean(dat$Age, na.rm=TRUE),
                                      sd(dat$Age, na.rm=TRUE)),
    "BMI (kg/m²), Median (IQR)" = sprintf("%.2f (%.2f–%.2f)",
                                          median(dat$BMI, na.rm=TRUE),
                                          quantile(dat$BMI, 0.25, na.rm=TRUE),
                                          quantile(dat$BMI, 0.75, na.rm=TRUE)),
    "Mode of Conception" = "",
    "Unassisted" = sprintf("%d (%.1f%%)",
                           sum(dat$ConceptionMode==1, na.rm=TRUE),
                           100*mean(dat$ConceptionMode==1, na.rm=TRUE)),
    "Medication and/or Insemination" = sprintf("%d (%.1f%%)",
                                               sum(dat$ConceptionMode==2, na.rm=TRUE),
                                               100*mean(dat$ConceptionMode==2, na.rm=TRUE)),
    "Assisted Reproductive Technology" = sprintf("%d (%.1f%%)",
                                                 sum(dat$ConceptionMode==3, na.rm=TRUE),
                                                 100*mean(dat$ConceptionMode==3, na.rm=TRUE)),
    "Race" = "",
    "White" = sprintf("%d (%.1f%%)", sum(dat$Race_Ethnicity==1, na.rm=TRUE),
                      100*mean(dat$Race_Ethnicity==1, na.rm=TRUE)),
    "Asian/PI" = sprintf("%d (%.1f%%)", sum(dat$Race_Ethnicity==2, na.rm=TRUE),
                         100*mean(dat$Race_Ethnicity==2, na.rm=TRUE)),
    "Black" = sprintf("%d (%.1f%%)", sum(dat$Race_Ethnicity==4, na.rm=TRUE),
                      100*mean(dat$Race_Ethnicity==4, na.rm=TRUE)),
    "Mixed or Other" = sprintf("%d (%.1f%%)", sum(dat$Race_Ethnicity==5, na.rm=TRUE),
                               100*mean(dat$Race_Ethnicity==5, na.rm=TRUE)),
    "Previous Live Birth" = "",
    "Yes" = sprintf("%d (%.1f%%)", sum(dat$Parous==1, na.rm=TRUE),
                    100*mean(dat$Parous==1, na.rm=TRUE)),
    "No" = sprintf("%d (%.1f%%)", sum(dat$Parous==0, na.rm=TRUE),
                   100*mean(dat$Parous==0, na.rm=TRUE))
  )
  tibble(Characteristic = names(out), !!group_name := unlist(out))
}

# --- Summaries ---
cases_summary <- summarize_group(cases, sprintf("Early Pregnancy Loss (N = %d)", nrow(cases)))
controls_summary <- summarize_group(controls, sprintf("Ongoing Pregnancy (N = %d)", nrow(controls)))

table1 <- full_join(cases_summary, controls_summary, by="Characteristic")

# --- Identify section headers ---
section_rows <- grep("Age|BMI|Mode of Conception|Race|Previous Live Birth", table1$Characteristic)

# --- Format table ---
gt_tbl <- table1 %>%
  gt() %>%
  tab_header(
    title = md("**Table 1. Participant demographics from the UCSF cohort (N = 69), stratified by case status.**")
  ) %>%
  cols_align(align="left", columns=everything()) %>%
  # Bold section headers
  tab_style(
    style = cell_text(weight="bold"),
    locations = cells_body(rows = section_rows)
  ) %>%
  # Table appearance — Times New Roman, no lines
  tab_options(
    table.font.size = 12,
    data_row.padding = px(4),
    column_labels.font.weight = "bold",
    table.font.names = "Times New Roman",
    table.border.top.width = px(0),
    table.border.bottom.width = px(0),
    table_body.hlines.width = px(0),
    table_body.vlines.width = px(0),
    heading.border.bottom.width = px(0),
    column_labels.border.top.width = px(0),
    column_labels.border.bottom.width = px(0),
    row_group.border.top.width = px(0),
    row_group.border.bottom.width = px(0)
  )

# --- Save ---
gtsave(gt_tbl, output_png, vwidth = 1000, vheight = 700)
cat("\n✅ Table 1 saved to:", output_png, "\n")
