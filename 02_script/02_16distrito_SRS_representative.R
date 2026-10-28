# ==============================================================================
# 專案：Zambézia 各縣 (Distrito) 抽樣點空間歸屬統計與互動地圖 (修正版)
# ==============================================================================

library(here)
library(sf)
library(dplyr)
library(mapview)

sf_use_s2(FALSE)

# ------------------------------------------------------------------------------
# 步驟 1：下載並篩選 Zambezia 縣級邊界 (Admin 2)
# ------------------------------------------------------------------------------
cat(">>> [1/5] 下載莫三比克標準縣級邊界 (GADM)...\n")

moz_adm2_raw <- geodata::gadm(country = "MOZ", level = 2, path = tempdir()) %>%
  st_as_sf()

# 使用無重音的 "Zambezia" 進行篩選
zambezia_distritos <- moz_adm2_raw %>%
  filter(NAME_1 == "Zambezia") %>%
  st_transform(4326) %>%
  st_make_valid() %>%
  select(
    district_gadm = NAME_2,
    geometry
  )

cat("  - Zambézia 轄下縣份總數 (Distritos):", nrow(zambezia_distritos), "\n")
cat("  - 縣名清單:\n   ", paste(sort(zambezia_distritos$district_gadm), collapse = ", "), "\n")


# ------------------------------------------------------------------------------
# 步驟 2：讀取 4,874 點位資料並執行空間交集
# ------------------------------------------------------------------------------
cat("\n>>> [2/5] 讀取 base_ac_processed.rds 並與縣邊界進行空間交集...\n")

base_ac_processed <- readRDS(here("01_data", "processed", "base_ac_processed.rds")) %>%
  st_transform(4326)

pts_with_district <- st_join(
  base_ac_processed,
  zambezia_distritos,
  join = st_intersects,
  left = TRUE
)

# 邊界外點位檢查 (若有邊界微小漂移，採用最近鄰自動校正)
na_count <- sum(is.na(pts_with_district$district_gadm))
if (na_count > 0) {
  cat("  - 偵測到", na_count, "個點位落在邊界外，啟動最近鄰校正...\n")
  na_indices <- which(is.na(pts_with_district$district_gadm))
  nearest_poly_idx <- st_nearest_feature(pts_with_district[na_indices, ], zambezia_distritos)
  pts_with_district$district_gadm[na_indices] <- zambezia_distritos$district_gadm[nearest_poly_idx]
}

cat("  - 點位歸位完成！未匹配點數 (NA):", sum(is.na(pts_with_district$district_gadm)), "\n")


# ------------------------------------------------------------------------------
# 步驟 3：樞紐分析統計各縣點位分佈
# ------------------------------------------------------------------------------
cat("\n>>> [3/5] 彙總各縣之抽樣點數 (4874 Frame, 115 SRS, 30 Sero)...\n")

district_summary_table <- pts_with_district %>%
  st_drop_geometry() %>%
  group_by(district_gadm) %>%
  summarize(
    Total_Frame_N4874 = n(),
    COMSA_SRS_N115    = sum(COMSA %in% c(1, 2)),
    COMSA_Sero_N30    = sum(COMSA == 2),
    SRS_Sampling_Rate = paste0(round(COMSA_SRS_N115 / Total_Frame_N4874 * 100, 2), "%"),
    Has_Sero_Sites    = if_else(COMSA_Sero_N30 > 0, "Yes", "No"),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Frame_N4874))

# 加上合計列
total_row <- data.frame(
  district_gadm     = "TOTAL (Zambézia)",
  Total_Frame_N4874 = sum(district_summary_table$Total_Frame_N4874),
  COMSA_SRS_N115    = sum(district_summary_table$COMSA_SRS_N115),
  COMSA_Sero_N30    = sum(district_summary_table$COMSA_Sero_N30),
  SRS_Sampling_Rate = paste0(round(sum(district_summary_table$COMSA_SRS_N115) / sum(district_summary_table$Total_Frame_N4874) * 100, 2), "%"),
  Has_Sero_Sites    = paste0(sum(district_summary_table$COMSA_Sero_N30 > 0), " / ", nrow(zambezia_distritos), " Districts")
)

final_summary_table <- bind_rows(district_summary_table, total_row)

cat("\n=== 各縣 (Distrito) 空間抽樣點位分佈統計表 ===\n")
print(as.data.frame(final_summary_table))

write.csv(
  final_summary_table,
  here("03_results", "table_district_sampling_summary.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 步驟 4：準備互動地圖視覺化物件
# ------------------------------------------------------------------------------
cat("\n>>> [4/5] 整理圖層屬性標籤...\n")

zambezia_distritos_vis <- zambezia_distritos %>%
  left_join(district_summary_table, by = "district_gadm") %>%
  mutate(
    Total_Frame_N4874 = coalesce(Total_Frame_N4874, 0L),
    COMSA_SRS_N115    = coalesce(COMSA_SRS_N115, 0L),
    COMSA_Sero_N30    = coalesce(COMSA_Sero_N30, 0L),
    SRS_Sampling_Rate = coalesce(SRS_Sampling_Rate, "0%"),
    Has_Sero_Sites    = coalesce(Has_Sero_Sites, "No"),
    tooltip_district  = paste0(
      "Distrito: ", district_gadm,
      " | Frame: ", Total_Frame_N4874,
      " | SRS: ", COMSA_SRS_N115,
      " | Sero: ", COMSA_Sero_N30
    )
  )

pts_vis <- pts_with_district %>%
  mutate(
    comsa_tier = case_when(
      COMSA == 2 ~ "COMSA + Serosurvey (N = 30)",
      COMSA == 1 ~ "COMSA Only (N = 85)",
      COMSA == 0 ~ "Baseline Sampling Frame (N = 4,759)",
      TRUE       ~ "Unclassified"
    ),
    comsa_tier = factor(
      comsa_tier,
      levels = c(
        "COMSA + Serosurvey (N = 30)",
        "COMSA Only (N = 85)",
        "Baseline Sampling Frame (N = 4,759)"
      )
    ),
    tooltip_pt = paste0(
      "Distrito: ", district_gadm,
      " | Tier: ", comsa_tier,
      if_else(!is.na(comsa_cluster_id), paste0(" | Cluster: ", comsa_cluster_id), "")
    )
  ) %>%
  arrange(COMSA)


# ------------------------------------------------------------------------------
# 步驟 5：渲染互動式地圖至 Viewer
# ------------------------------------------------------------------------------
cat("\n>>> [5/5] 渲染互動地圖至 RStudio Viewer...\n")

mapviewOptions(
  basemaps = c("OpenStreetMap", "Esri.WorldImagery", "CartoDB.Positron"),
  fgb = TRUE
)

tier_colors <- c(
  "COMSA + Serosurvey (N = 30)"         = "#00E676",
  "COMSA Only (N = 85)"                = "#FF9100",
  "Baseline Sampling Frame (N = 4,759)"= "#2979FF"
)

# 縣界圖層 (原生 Popup)
map_district_layer <- mapview(
  zambezia_distritos_vis,
  color = "#222222",
  lwd = 1.6,
  alpha.regions = 0.08,
  label = zambezia_distritos_vis$tooltip_district,
  layer.name = "Zambézia Distritos (Admin 2)",
  popup = TRUE
)

# 點位圖層 (原生 Popup)
map_points_layer <- mapview(
  pts_vis,
  zcol = "comsa_tier",
  col.regions = tier_colors,
  color = "black",
  lwd = 0.5,
  cex = 3.5,
  alpha.regions = 0.85,
  label = pts_vis$tooltip_pt,
  layer.name = "Sampling Tiers (N = 4,874)",
  popup = TRUE
)

map_district_view <- map_district_layer + map_points_layer

sf_use_s2(TRUE)
print(map_district_view)



# Inspect Data

# ------------------------------------------------------------------------------
# 診斷 1：檢查 base_ac_processed 原始標籤到底有幾個
# ------------------------------------------------------------------------------
cat("=== 檢查 base_ac_processed 本身的分層計數 ===\n")
table(base_ac_processed$COMSA, useNA = "always")

# 檢查指派到的唯一 cluster 數量
assigned_clusters <- unique(na.omit(base_ac_processed$comsa_cluster_id))
cat("成功指派的 COMSA Cluster 總數:", length(assigned_clusters), "(預期應為 115)\n")

# ------------------------------------------------------------------------------
# 診斷 2：找出到底是哪 7 個 Cluster 沒有被選中
# ------------------------------------------------------------------------------
all_comsa_clusters <- comsa_dissolved$cluster
missing_clusters   <- setdiff(all_comsa_clusters, assigned_clusters)

cat("\n未成功在 4,874 點位中找到代表點的 Cluster ID:\n")
print(missing_clusters)

# 檢查這 7 個缺失 Cluster 的 sero 屬性 (是否包含那個丟失的 sero == 1)
missing_info <- comsa_dissolved %>%
  st_drop_geometry() %>%
  filter(cluster %in% missing_clusters) %>%
  select(cluster, distrito, sero, any_of("seroprev"))

cat("\n缺失 Cluster 的詳細屬性（檢查 sero 是否為 1）：\n")
print(missing_info)

# ------------------------------------------------------------------------------
# 診斷 3：檢查這 7 個多邊形內部是否真的「0 個點」
# ------------------------------------------------------------------------------
# 抽出缺失的 7 個多邊形
missing_polys <- comsa_dissolved %>% filter(cluster %in% missing_clusters)

# 重新計算落在這 7 個多邊形內的點位數
pts_in_missing <- st_join(base_ac_sf, missing_polys, join = st_intersects, left = FALSE)
cat("\n落入這 7 個多邊形內部的 Base AC 點位數:", nrow(pts_in_missing), "\n")