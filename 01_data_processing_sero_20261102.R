# ==============================================================================
# 專案：莫三比克 Zambézia 省空間抽樣資料集整併與 15 種抗原血清盛行率轉錄 (批次完整版)
# 產出檔案：01_data/processed/base_ac_processed.rds
# ==============================================================================

library(here)
library(sf)
library(dplyr)
library(tidyr)
library(mapview)
library(leafpop)

# 全域幾何設定：關閉 s2 引擎，改採 GEOS 平面幾何避免自交報錯
sf_use_s2(FALSE)

# ------------------------------------------------------------------------------
# 步驟 1：資料前處理與拓撲融合 (COMSA 115 Clusters & Centroids)
# ------------------------------------------------------------------------------
cat(">>> 執行步驟 1：COMSA EA 多邊形修復、拓撲融合與幾何中心摘要...\n")

comsa_raw <- st_read(
  here("01_data", "raw", "comsa.zambezia_with_clusters", "comsa.zambezia_with_clusters.shp"),
  quiet = TRUE
) %>% st_transform(4326)

comsa_dissolved <- comsa_raw %>%
  filter(cluster != 999) %>%
  st_make_valid() %>%
  group_by(cluster) %>%
  summarize(
    sero = max(sero, na.rm = TRUE),
    distrito = first(distrito),
    geometry = st_union(st_combine(geometry)),
    .groups = "drop"
  ) %>%
  st_make_valid()

summarized_comsa_sf <- comsa_dissolved %>%
  st_centroid()

cat("  - COMSA 融合後 Cluster 總數:", nrow(comsa_dissolved), "\n")
cat("  - 幾何中心點數量:", nrow(summarized_comsa_sf), "\n")


# ------------------------------------------------------------------------------
# 步驟 2：Mozambique 血清資料家戶幾何與「15 種抗原」血清陽性率批次計算
# ------------------------------------------------------------------------------
cat("\n>>> 執行步驟 2：批次計算 Moz Age 30 個群聚之幾何中心與 15 種抗原 seroprev...\n")

moz_age <- readRDS(here("01_data", "raw", "moz_age.rds"))

# 2.1 幾何中心：以 hhid 去重後計算經緯度平均
moz_cluster_coords <- moz_age %>%
  group_by(csa_code, hhid) %>%
  summarize(
    lon = first(lon),
    lat = first(lat),
    .groups = "drop"
  ) %>%
  group_by(csa_code) %>%
  summarize(
    mean_lon = mean(lon, na.rm = TRUE),
    mean_lat = mean(lat, na.rm = TRUE),
    n_households = n(),
    .groups = "drop"
  )

# 2.2 自動辨識並批次計算全部 15 種 pos_* 抗原盛行率
pos_cols <- grep("^pos_", names(moz_age), value = TRUE)
cat("  - 偵測到之抗原清單 (共", length(pos_cols), "種):\n   ", paste(pos_cols, collapse = ", "), "\n")

moz_cluster_sero <- moz_age %>%
  group_by(csa_code) %>%
  summarize(
    n_individuals = n(),
    # 批次將 pos_* 轉為 seroprev_*
    across(
      all_of(pos_cols),
      ~ mean(.x, na.rm = TRUE),
      .names = "seroprev_{sub('^pos_', '', .col)}"
    ),
    .groups = "drop"
  )

# 2.3 整併為 30 個群聚的空間資料物件
summarized_moz_age_sf <- moz_cluster_coords %>%
  left_join(moz_cluster_sero, by = "csa_code") %>%
  st_as_sf(coords = c("mean_lon", "mean_lat"), crs = 4326, remove = FALSE) %>%
  mutate(
    tooltip_moz = paste0("CSA: ", csa_code, " | Indiv: ", n_individuals, " (HH: ", n_households, ")")
  )

cat("  - Moz Age 彙總群聚數:", nrow(summarized_moz_age_sf), "\n")


# ------------------------------------------------------------------------------
# 步驟 3：互動式衛星地圖檢驗
# ------------------------------------------------------------------------------
cat("\n>>> 執行步驟 3：產出互動式地圖以檢核 30 個 Centroids 分佈合理性...\n")

comsa_ea_vis <- comsa_raw %>%
  st_transform(4326) %>%
  mutate(
    sero_status = case_when(
      sero == 1 ~ "Selected (Serosurvey)",
      sero == 0 ~ "Not Selected",
      TRUE      ~ "Unknown"
    ),
    sero_status = factor(sero_status, levels = c("Selected (Serosurvey)", "Not Selected")),
    tooltip_ea = paste0("EA Cluster: ", cluster, " | Status: ", sero_status)
  )

mapviewOptions(
  basemaps = c("Esri.WorldImagery", "OpenStreetMap"),
  fgb = TRUE
)

map_check <- mapview(
  comsa_ea_vis,
  zcol = "sero_status",
  col.regions = c("Selected (Serosurvey)" = "#FF3333", "Not Selected" = "#33CCFF"),
  color = "white",
  lwd = 1.0,
  alpha.regions = 0.35,
  label = comsa_ea_vis$tooltip_ea,
  layer.name = "COMSA EAs (341 Polygons)"
) +
  mapview(
    summarized_moz_age_sf,
    col.regions = "#FFFF00",
    color = "black",
    lwd = 1.2,
    cex = 6,
    alpha.regions = 0.95,
    label = summarized_moz_age_sf$tooltip_moz,
    layer.name = "Moz 30 Centroids"
  )

print(map_check)
cat("  - 互動地圖已成功輸出至 Viewer。\n")


# ------------------------------------------------------------------------------
# 步驟 4：大圓距離配對與 15 種抗原 seroprev 轉錄 (Moz -> COMSA)
# ------------------------------------------------------------------------------
cat("\n>>> 執行步驟 4：透過最近鄰大圓距離配對轉錄 15 種抗原盛行率...\n")

sf_use_s2(TRUE)
dist_matrix <- st_distance(summarized_moz_age_sf, summarized_comsa_sf)
sf_use_s2(FALSE)

nearest_comsa_indices <- apply(dist_matrix, 1, which.min)
nearest_distances_km  <- apply(dist_matrix, 1, min) / 1000

# 建立包含全部 seroprev_* 的配對資料表
seroprev_cols <- grep("^seroprev_", names(summarized_moz_age_sf), value = TRUE)

matching_df <- data.frame(
  csa_code = summarized_moz_age_sf$csa_code,
  matched_comsa_cluster = summarized_comsa_sf$cluster[nearest_comsa_indices],
  dist_km = nearest_distances_km,
  n_individuals = summarized_moz_age_sf$n_individuals
) %>%
  bind_cols(st_drop_geometry(summarized_moz_age_sf)[, seroprev_cols])

cat("  - 配對最大距離 (Max distance):", round(max(matching_df$dist_km), 2), "km\n")
cat("  - 配對中位數距離 (Median distance):", round(median(matching_df$dist_km), 2), "km\n")

if (length(unique(matching_df$matched_comsa_cluster)) != 30) {
  warning("警告：存在多個 Moz Cluster 配對到同一個 COMSA Cluster！")
} else {
  cat("  - 驗證成功：30 個 Moz Clusters 與 COMSA 達成嚴格 1 對 1 配對。\n")
}

# 轉錄進 COMSA 多邊形與 Centroid 物件
cols_to_merge <- matching_df %>% select(matched_comsa_cluster, n_individuals, all_of(seroprev_cols))

comsa_dissolved <- comsa_dissolved %>%
  left_join(cols_to_merge, by = c("cluster" = "matched_comsa_cluster"))

summarized_comsa_sf <- summarized_comsa_sf %>%
  left_join(cols_to_merge, by = c("cluster" = "matched_comsa_cluster"))


# ------------------------------------------------------------------------------
# 步驟 5：Base AC 4,874 點位空間包含檢驗、代表點裁決與完整抗原轉錄
# ------------------------------------------------------------------------------
cat("\n>>> 執行步驟 5：Base AC 4,874 點位空間裁決與全部抗原欄位掛載...\n")

base_ac_sf <- st_read(
  here("01_data", "raw", "ZBAC_Points_Final", "04_Base_AC_Final.shp"),
  quiet = TRUE
) %>% st_transform(4326)

base_ac_sf$temp_row_id <- seq_len(nrow(base_ac_sf))

# 空間包含檢驗
poly_sub <- comsa_dissolved %>%
  select(comsa_cluster = cluster, poly_sero = sero, n_individuals, all_of(seroprev_cols))

pts_in_poly <- st_join(
  base_ac_sf,
  poly_sub,
  join = st_intersects,
  left = FALSE
)

# 最近中心點距離計算與唯一裁決
pts_with_dist <- pts_in_poly %>%
  left_join(
    st_drop_geometry(summarized_comsa_sf) %>% select(cluster),
    by = c("comsa_cluster" = "cluster")
  )

candidate_distances <- numeric(nrow(pts_with_dist))
for (i in seq_len(nrow(pts_with_dist))) {
  target_cluster <- pts_with_dist$comsa_cluster[i]
  centroid_geom <- st_geometry(summarized_comsa_sf[summarized_comsa_sf$cluster == target_cluster, ])
  candidate_distances[i] <- as.numeric(st_distance(st_geometry(pts_with_dist[i, ]), centroid_geom))
}
pts_with_dist$dist_to_centroid <- candidate_distances

primary_selected_pts <- pts_with_dist %>%
  group_by(comsa_cluster) %>%
  slice_min(order_by = dist_to_centroid, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("  - 經最近中心點裁決後，確認指派的 COMSA 代表點數量:", nrow(primary_selected_pts), "\n")

# 整併回全體 4,874 點位
base_ac_processed <- base_ac_sf %>%
  left_join(
    st_drop_geometry(primary_selected_pts) %>%
      select(temp_row_id, comsa_cluster, poly_sero, n_individuals, all_of(seroprev_cols)),
    by = "temp_row_id"
  ) %>%
  mutate(
    COMSA = case_when(
      !is.na(comsa_cluster) & poly_sero == 1 ~ 2,
      !is.na(comsa_cluster) & poly_sero == 0 ~ 1,
      TRUE                                   ~ 0
    ),
    comsa_cluster_id = comsa_cluster,
    n_sample_size = if_else(COMSA == 2, n_individuals, NA_integer_)
  ) %>%
  # 僅在 COMSA == 2 的 30 個點保留血清盛行率數值，其餘為 NA
  mutate(
    across(
      all_of(seroprev_cols),
      ~ if_else(COMSA == 2, .x, NA_real_)
    )
  ) %>%
  select(-temp_row_id, -comsa_cluster, -poly_sero, -n_individuals)


# ------------------------------------------------------------------------------
# 步驟 6：完整性驗證與 RDS 持久化儲存
# ------------------------------------------------------------------------------
cat("\n>>> 執行步驟 6：資料完整性驗證與 RDS 匯出...\n")

sf_use_s2(TRUE)

cat("  - 總筆數確認 (應為 4874):", nrow(base_ac_processed), "\n")
cat("  - COMSA 分類點位數分佈：\n")
print(table(COMSA = base_ac_processed$COMSA, useNA = "always"))

cat("\n  - 各抗原 seroprev 具備數值的筆數 (應精確為 30 筆)：\n")
non_na_counts <- sapply(base_ac_processed %>% st_drop_geometry() %>% select(all_of(seroprev_cols)), function(x) sum(!is.na(x)))
print(non_na_counts)

# 儲存
output_dir <- here("01_data", "processed")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_path <- here("01_data", "processed", "base_ac_processed.rds")
saveRDS(base_ac_processed, file = output_path)

cat("\n  ✓ 包含 15 種抗原 seroprev 的主資料集已安全儲存至：\n   ", output_path, "\n")