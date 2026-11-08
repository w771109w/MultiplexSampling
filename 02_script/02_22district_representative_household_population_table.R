# ==============================================================================
# 專案：Zambézia 抽樣母體與子集之人口與家戶分層統計報表 (含 1,292 實際血清樣本數)
# 輸入：01_data/processed/base_ac_processed.rds
# 輸出：
#   - 03_results/table1_district_population_coverage.csv
#   - 03_results/table2_district_household_tiers.csv
# ==============================================================================

library(here)
library(sf)
library(dplyr)
library(tidyr)
library(scales)

# 確保輸出目錄存在
dir.create(here("03_results"), recursive = TRUE, showWarnings = FALSE)

# 關閉 s2 平面運算以確保多邊形融合拓撲穩定
sf_use_s2(FALSE)

# ------------------------------------------------------------------------------
# 步驟 1：建構 22 縣邊界並透過空間交集統一行政區定義
# ------------------------------------------------------------------------------
cat(">>> [1/3] 解析 GADM Level 3 重組現行 22 縣邊界，並與點位執行空間交集...\n")

if (!requireNamespace("geodata", quietly = TRUE)) install.packages("geodata")
moz_adm3_raw <- geodata::gadm(country = "MOZ", level = 3, path = tempdir()) %>%
  st_as_sf()

# 利用 GADM Level 3 (Postos) 精確重組為改制後的現行 22 縣
zambezia_distritos <- moz_adm3_raw %>%
  filter(NAME_1 == "Zambezia") %>%
  mutate(
    district_22 = case_when(
      NAME_3 == "Derre"                                       ~ "Derre",
      NAME_3 == "Luabo"                                       ~ "Luabo",
      NAME_3 %in% c("Mocubela", "Bajone")                     ~ "Mocubela",
      NAME_3 %in% c("Molumbo", "Corromana")                   ~ "Molumbo",
      grepl("^Mulevala", NAME_3, ignore.case = TRUE)           ~ "Mulevala",
      NAME_2 == "Nicoadala" & grepl("Quelimane", NAME_3)         ~ "Quelimane",
      TRUE                                                    ~ NAME_2
    )
  ) %>%
  st_make_valid() %>%
  group_by(district_22) %>%
  summarize(
    geometry = st_union(st_combine(geometry)),
    .groups = "drop"
  ) %>%
  st_make_valid() %>%
  st_transform(4326) %>%
  rename(district_gadm = district_22)

# 讀取 base_ac 空間物件
base_ac_raw <- readRDS(here("01_data", "processed", "base_ac_processed.rds")) %>%
  st_transform(4326)

# 若原本帶有舊的 district_gadm 欄位，先剔除避免干擾
if ("district_gadm" %in% names(base_ac_raw)) {
  base_ac_raw <- base_ac_raw %>% select(-district_gadm)
}

# 執行空間交集，以多邊形邊界賦予真正的 22 縣標籤
pts_with_district <- st_join(
  base_ac_raw,
  zambezia_distritos,
  join = st_intersects,
  left = TRUE
)

# 邊界外點位檢查 (若有沿海微小漂移，以最近鄰校正)
na_indices <- which(is.na(pts_with_district$district_gadm))
if (length(na_indices) > 0) {
  nearest_poly_idx <- st_nearest_feature(pts_with_district[na_indices, ], zambezia_distritos)
  pts_with_district$district_gadm[na_indices] <- zambezia_distritos$district_gadm[nearest_poly_idx]
}

# 移除空間幾何進入統計分析
base_ac <- pts_with_district %>% st_drop_geometry()

# 自動對齊家戶數、人口、實際樣本數欄位
hh_col <- if ("N_AF_AC" %in% names(base_ac)) "N_AF_AC" else "AF"
pop_col <- if ("N_THM_AC" %in% names(base_ac)) "N_THM_AC" else "THM"
sample_col <- if ("n_sample_size" %in% names(base_ac)) "n_sample_size" else "n_individuals"
district_col <- "district_gadm"

cat("  - 使用行政區欄位: district_gadm (已對齊 22 縣空間邊界)\n")
cat("  - 使用人口欄位:", pop_col, "\n")
cat("  - 使用家戶欄位:", hh_col, "\n")
cat("  - 使用實際樣本數欄位:", sample_col, "\n")


# ------------------------------------------------------------------------------
# 步驟 2：產出表格 1 —— 各縣人口數 (Persons)、Sero / SRS 覆蓋佔比及實際血清樣本數
# ------------------------------------------------------------------------------
cat("\n>>> [2/3] 建立表格 1：各縣總人口、抽樣覆蓋佔比與實際血清樣本數 (N=1,292)...\n")

table1_district <- base_ac %>%
  group_by(district = .data[[district_col]]) %>%
  summarize(
    total_pop       = sum(.data[[pop_col]], na.rm = TRUE),
    srs_pop         = sum(if_else(COMSA %in% c(1, 2), .data[[pop_col]], 0), na.rm = TRUE),
    sero_pop        = sum(if_else(COMSA == 2, .data[[pop_col]], 0), na.rm = TRUE),
    sample_pop      = sum(if_else(COMSA == 2, as.numeric(.data[[sample_col]]), 0), na.rm = TRUE),
    .groups = "drop"
  )

# 確保 22 縣完整列出（即使該縣母體或抽樣為 0）
all_districts_df <- tibble(district = zambezia_distritos$district_gadm)
table1_district <- all_districts_df %>%
  left_join(table1_district, by = "district") %>%
  mutate(
    total_pop   = coalesce(total_pop, 0),
    srs_pop     = coalesce(srs_pop, 0),
    sero_pop    = coalesce(sero_pop, 0),
    sample_pop  = coalesce(sample_pop, 0),
    srs_pop_pct    = if_else(total_pop > 0, srs_pop / total_pop, 0),
    sero_pop_pct   = if_else(total_pop > 0, sero_pop / total_pop, 0),
    sample_pop_pct = if_else(total_pop > 0, sample_pop / total_pop, 0)
  ) %>%
  arrange(desc(total_pop))

# 計算全省加總總計列 (Total Row)
table1_total <- tibble(
  district        = "Total (Zambézia)",
  total_pop       = sum(table1_district$total_pop),
  srs_pop         = sum(table1_district$srs_pop),
  sero_pop        = sum(table1_district$sero_pop),
  sample_pop      = sum(table1_district$sample_pop),
  srs_pop_pct     = srs_pop / total_pop,
  sero_pop_pct    = sero_pop / total_pop,
  sample_pop_pct  = sample_pop / total_pop
)

# 合併並格式化
table1_final <- bind_rows(table1_district, table1_total) %>%
  mutate(
    `Total Pop`       = comma(total_pop),
    `COMSA SRS Pop`   = comma(srs_pop),
    `SRS Pop (%)`     = percent(srs_pop_pct, accuracy = 0.1),
    `Sero Pop`        = comma(sero_pop),
    `Sero Pop (%)`    = percent(sero_pop_pct, accuracy = 0.1),
    `Sero Sample (N)` = comma(sample_pop),
    `Sero Sample (%)` = percent(sample_pop_pct, accuracy = 0.001)
  ) %>%
  select(
    District = district,
    `Total Pop`,
    `COMSA SRS Pop`,
    `SRS Pop (%)`,
    `Sero Pop`,
    `Sero Pop (%)`,
    `Sero Sample (N)`,
    `Sero Sample (%)`
  )

# 輸出 CSV 與終端檢視
write.csv(table1_final, here("03_results", "table1_district_population_coverage.csv"), row.names = FALSE)
cat("--- 表格 1 預覽 (更新欄位名稱) ---\n")
print(as.data.frame(table1_final))


# ------------------------------------------------------------------------------
# 步驟 3：產出表格 2 —— 各縣家戶數 (Households) 於三階層之分佈與群聚數
# ------------------------------------------------------------------------------
cat("\n>>> [3/3] 建立表格 2：各縣家戶數跨抽樣階層比較...\n")

table2_district <- base_ac %>%
  group_by(district = .data[[district_col]]) %>%
  summarize(
    # 全母體層
    frame_clusters  = n(),
    frame_hh        = sum(.data[[hh_col]], na.rm = TRUE),
    
    # COMSA SRS 層
    srs_clusters    = sum(COMSA %in% c(1, 2)),
    srs_hh          = sum(if_else(COMSA %in% c(1, 2), .data[[hh_col]], 0), na.rm = TRUE),
    
    # COMSA Sero 層
    sero_clusters   = sum(COMSA == 2),
    sero_hh         = sum(if_else(COMSA == 2, .data[[hh_col]], 0), na.rm = TRUE),
    .groups = "drop"
  )

# 補齊 22 縣完整列
table2_district <- all_districts_df %>%
  left_join(table2_district, by = "district") %>%
  mutate(
    frame_clusters = coalesce(frame_clusters, 0L),
    frame_hh       = coalesce(frame_hh, 0),
    srs_clusters   = coalesce(srs_clusters, 0L),
    srs_hh         = coalesce(srs_hh, 0),
    sero_clusters  = coalesce(sero_clusters, 0L),
    sero_hh        = coalesce(sero_hh, 0),
    srs_hh_pct     = if_else(frame_hh > 0, srs_hh / frame_hh, 0),
    sero_hh_pct    = if_else(frame_hh > 0, sero_hh / frame_hh, 0)
  ) %>%
  arrange(desc(frame_hh))

# 計算全省總計列
table2_total <- tibble(
  district        = "Total (Zambézia)",
  frame_clusters  = sum(table2_district$frame_clusters),
  frame_hh        = sum(table2_district$frame_hh),
  srs_clusters    = sum(table2_district$srs_clusters),
  srs_hh          = sum(table2_district$srs_hh),
  sero_clusters   = sum(table2_district$sero_clusters),
  sero_hh         = sum(table2_district$sero_hh),
  srs_hh_pct      = srs_hh / frame_hh,
  sero_hh_pct     = sero_hh / frame_hh
)

# 格式化輸出
table2_final <- bind_rows(table2_district, table2_total) %>%
  mutate(
    `Frame Clusters`  = comma(frame_clusters),
    `Frame HH`        = comma(frame_hh),
    `SRS Clusters`    = comma(srs_clusters),
    `SRS HH`          = comma(srs_hh),
    `SRS HH (%)`      = percent(srs_hh_pct, accuracy = 0.1),
    `Sero Clusters`   = comma(sero_clusters),
    `Sero HH`         = comma(sero_hh),
    `Sero HH (%)`     = percent(sero_hh_pct, accuracy = 0.1)
  ) %>%
  select(
    District = district,
    `Frame Clusters`, `Frame HH`,
    `SRS Clusters`, `SRS HH`, `SRS HH (%)`,
    `Sero Clusters`, `Sero HH`, `Sero HH (%)`
  )

# 輸出 CSV 與終端檢視
write.csv(table2_final, here("03_results", "table2_district_household_tiers.csv"), row.names = FALSE)
cat("\n--- 表格 2 預覽 (家戶數與群聚數跨階層分佈) ---\n")
print(as.data.frame(table2_final))

sf_use_s2(TRUE)
cat("\n✓ 兩張表格已順利建立並存入: 03_results/\n")