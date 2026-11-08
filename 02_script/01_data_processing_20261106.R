# ==============================================================================
# 專案：莫三比克 Zambézia 省空間抽樣資料集整併與 15 種抗原血清盛行率轉錄 (兩階段配對完整版)
# 產出檔案：01_data/processed/base_ac_processed.rds
# 特點：
#   1. 階段一：Base AC 點位優先以多邊形空間包含 (st_intersects) 指派離 Centroid 最近點
#   2. 階段二：針對多邊形內無點位的 Cluster，啟動 Centroid 最近鄰搜尋進行救援補齊
#   3. 確保 30 個 Sero Clusters 與全數 SRS Clusters 完整被指定，無任何缺漏
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
# 步驟 1：資料前處理與拓撲融合 (COMSA Clusters & Centroids)
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
cat("  - 偵測到之抗原清單 (共", length(pos_cols), "種):\n    ", paste(pos_cols, collapse = ", "), "\n")

moz_cluster_sero <- moz_age %>%
  group_by(csa_code) %>%
  summarize(
    n_individuals = n(),
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
  layer.name = "COMSA EAs"
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
# 步驟 5：兩階段代表點裁決 (Polygon 內優先，Centroid 最近鄰救援)
# ------------------------------------------------------------------------------
cat("\n>>> 執行步驟 5：兩階段代表點裁決 (Polygon 內優先，Centroid 最近鄰救援)...\n")

base_ac_sf <- st_read(
  here("01_data", "raw", "ZBAC_Points_Final", "04_Base_AC_Final.shp"),
  quiet = TRUE
) %>% st_transform(4326)

base_ac_sf$temp_row_id <- seq_len(nrow(base_ac_sf))

# 準備 COMSA 屬性對照表
poly_sub <- comsa_dissolved %>%
  select(comsa_cluster = cluster, poly_sero = sero, n_individuals, all_of(seroprev_cols))

# ----------------- 階段 1：多邊形空間包含 (st_intersects) -----------------
pts_in_poly <- st_join(
  base_ac_sf,
  poly_sub,
  join = st_intersects,
  left = FALSE
)

# 計算多邊形內候選點到其 Cluster Centroid 的大圓距離
sf_use_s2(TRUE)
p1_dists <- numeric(nrow(pts_in_poly))
for (i in seq_len(nrow(pts_in_poly))) {
  target_cluster <- pts_in_poly$comsa_cluster[i]
  c_geom <- st_geometry(summarized_comsa_sf[summarized_comsa_sf$cluster == target_cluster, ])
  p1_dists[i] <- as.numeric(st_distance(st_geometry(pts_in_poly[i, ]), c_geom))
}
pts_in_poly$dist_to_centroid <- p1_dists
sf_use_s2(FALSE)

# 每個 cluster 在其多邊形內取距離 Centroid 最近的點
primary_p1 <- pts_in_poly %>%
  group_by(comsa_cluster) %>%
  slice_min(order_by = dist_to_centroid, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("  - 階段 1 (多邊形內命中指派):", nrow(primary_p1), "個 clusters\n")

# ----------------- 階段 2：針對遺漏 Clusters 執行 Centroid 救援 -----------------
all_comsa_clusters <- summarized_comsa_sf$cluster
missing_clusters <- setdiff(all_comsa_clusters, primary_p1$comsa_cluster)

if (length(missing_clusters) > 0) {
  cat("  - 階段 2 (啟動救援): 共有", length(missing_clusters), "個 clusters 多邊形內無點位，進行 Centroid 最近鄰搜尋...\n")
  
  # 排除階段 1 已選定的點位，避免同一個 Base AC 點重複指派
  available_base_sf <- base_ac_sf %>%
    filter(!temp_row_id %in% primary_p1$temp_row_id)
  
  fallback_list <- list()
  
  sf_use_s2(TRUE)
  for (cl_id in missing_clusters) {
    cl_centroid <- summarized_comsa_sf %>% filter(cluster == cl_id)
    
    # 計算到所有可用候選點的距離
    dist_vec <- as.numeric(st_distance(cl_centroid, available_base_sf))
    min_idx <- which.min(dist_vec)
    chosen_pt <- available_base_sf[min_idx, ]
    
    cl_attr <- st_drop_geometry(comsa_dissolved %>% filter(cluster == cl_id))
    
    fallback_row <- tibble(
      temp_row_id      = chosen_pt$temp_row_id,
      comsa_cluster    = cl_id,
      poly_sero        = cl_attr$sero,
      n_individuals    = cl_attr$n_individuals,
      dist_to_centroid = dist_vec[min_idx]
    ) %>%
      bind_cols(cl_attr %>% select(all_of(seroprev_cols)))
    
    fallback_list[[as.character(cl_id)]] <- fallback_row
    
    # 該點指派後自可用池中剔除
    available_base_sf <- available_base_sf[-min_idx, ]
  }
  sf_use_s2(FALSE)
  
  primary_p2 <- bind_rows(fallback_list)
  cat("  - 階段 2 救援配對完成，平均搜尋距離:", round(mean(primary_p2$dist_to_centroid) / 1000, 2), "km\n")
  
  # 合併階段 1 與階段 2 裁決結果
  final_selected_pts <- bind_rows(
    primary_p1 %>% st_drop_geometry() %>% select(temp_row_id, comsa_cluster, poly_sero, n_individuals, all_of(seroprev_cols)),
    primary_p2 %>% select(temp_row_id, comsa_cluster, poly_sero, n_individuals, all_of(seroprev_cols))
  )
} else {
  final_selected_pts <- primary_p1 %>% 
    st_drop_geometry() %>% 
    select(temp_row_id, comsa_cluster, poly_sero, n_individuals, all_of(seroprev_cols))
}

cat("  - 最終裁決代表點總數:", nrow(final_selected_pts), "\n")
cat("  - 其中 Sero 代表點數量 (應精確為 30):", sum(final_selected_pts$poly_sero == 1), "\n")

# ----------------- 整併回全體 4,874 點位 -----------------
base_ac_processed <- base_ac_sf %>%
  left_join(final_selected_pts, by = "temp_row_id") %>%
  mutate(
    COMSA = case_when(
      !is.na(comsa_cluster) & poly_sero == 1 ~ 2,
      !is.na(comsa_cluster) & poly_sero == 0 ~ 1,
      TRUE                                  ~ 0
    ),
    comsa_cluster_id = comsa_cluster,
    n_sample_size = if_else(COMSA == 2, n_individuals, NA_integer_)
  ) %>%
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
cat("  - COMSA 分類點位數分佈 (應確認 COMSA==2 為 30 筆)：\n")
print(table(COMSA = base_ac_processed$COMSA, useNA = "always"))

cat("\n  - 各抗原 seroprev 具備數值的筆數 (應精確為 30 筆)：\n")
non_na_counts <- sapply(
  base_ac_processed %>% st_drop_geometry() %>% select(all_of(seroprev_cols)), 
  function(x) sum(!is.na(x))
)
print(non_na_counts)

# 儲存輸出
output_dir <- here("01_data", "processed")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_path <- here("01_data", "processed", "base_ac_processed.rds")
saveRDS(base_ac_processed, file = output_path)

cat("\n  ✓ 包含 15 種抗原 seroprev 的主資料集已安全儲存至：\n    ", output_path, "\n")