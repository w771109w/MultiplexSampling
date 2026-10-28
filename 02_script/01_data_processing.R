# ==============================================================================
# 專案：莫三比克 Zambézia 省空間抽樣資料集整併與血清盛行率轉錄 (含步驟 3 互動地圖檢驗)
# 產出檔案：01_data/processed/base_ac_processed.rds
# ==============================================================================

library(here)
library(sf)
library(dplyr)
library(tidyr)
library(mapview)
library(leafpop)

# 全域幾何設定：關閉 s2 球面幾何引擎，改採 GEOS 平面幾何引擎以杜絕多邊形自交報錯
sf_use_s2(FALSE)

# ------------------------------------------------------------------------------
# 步驟 1：資料前處理與拓撲融合 (COMSA 115 Clusters & Centroids)
# ------------------------------------------------------------------------------
cat(">>> 執行步驟 1：COMSA EA 多邊形修復、拓撲融合與幾何中心摘要...\n")

# 讀取 COMSA 原始 EA 多邊形 (341 個多邊形)
comsa_raw <- st_read(
  here("01_data", "raw", "comsa.zambezia_with_clusters", "comsa.zambezia_with_clusters.shp"),
  quiet = TRUE
) %>% st_transform(4326)

# 排除 999 雜訊，修復拓撲自相交，依 cluster 融合成 115 個複合多邊形
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

# 計算 115 個 Cluster 的幾何中心點 (Centroids)
summarized_comsa_sf <- comsa_dissolved %>%
  st_centroid()

cat("  - COMSA 融合後 Cluster 總數:", nrow(comsa_dissolved), "\n")
cat("  - 幾何中心點數量:", nrow(summarized_comsa_sf), "\n")


# ------------------------------------------------------------------------------
# 步驟 2：Mozambique 血清資料家戶幾何與血清陽性率彙總 (30 Clusters)
# ------------------------------------------------------------------------------
cat("\n>>> 執行步驟 2：計算 Moz Age 30 個群聚之幾何中心與 seroprev...\n")

moz_age <- readRDS(here("01_data", "raw", "moz_age.rds"))

# 2.1 幾何中心：先以家戶 (hhid) 為單位去重，再算各 CSA 群聚的經緯度算術平均
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

# 2.2 血清盛行率：在個人層級 (Individual level) 計算抗體陽性率
target_pos_col <- grep("^pos", names(moz_age), value = TRUE)[1]
if (is.na(target_pos_col)) target_pos_col <- "pos" # 防呆欄位

moz_cluster_sero <- moz_age %>%
  group_by(csa_code) %>%
  summarize(
    n_individuals = n(),
    seroprev = mean(.data[[target_pos_col]], na.rm = TRUE),
    .groups = "drop"
  )

# 2.3 整併為 30 個群聚的空間資料物件
summarized_moz_age_sf <- moz_cluster_coords %>%
  left_join(moz_cluster_sero, by = "csa_code") %>%
  st_as_sf(coords = c("mean_lon", "mean_lat"), crs = 4326, remove = FALSE) %>%
  mutate(
    tooltip_moz = paste0("CSA: ", csa_code, " | Seroprev: ", round(seroprev, 3), " (HH: ", n_households, ")")
  )

cat("  - Moz Age 彙總群聚數:", nrow(summarized_moz_age_sf), "\n")


# ------------------------------------------------------------------------------
# 步驟 3：互動式衛星地圖檢驗 (檢查 30 個 Moz Centroids 是否精準落於 sero==1 EAs)
# ------------------------------------------------------------------------------
cat("\n>>> 執行步驟 3：產出互動式地圖以檢核 30 個 Centroids 分佈合理性...\n")

# 準備 COMSA 341 個 EA 多邊形以利檢視微觀邊界
comsa_ea_vis <- comsa_raw %>%
  st_transform(4326) %>%
  mutate(
    sero_status = case_when(
      sero == 1 ~ "Selected (Serosurvey)",
      sero == 0 ~ "Not Selected",
      TRUE      ~ "Unknown"
    ),
    sero_status = factor(sero_status, levels = c("Selected (Serosurvey)", "Not Selected")),
    tooltip_ea = paste0("EA Cluster: ", cluster, " | Status: ", sero_status, " (Row: ", row, ")")
  )

# 設定免 API Key 圖磚：優先載入真實高解析衛星空照圖
mapviewOptions(
  basemaps = c("Esri.WorldImagery", "OpenStreetMap"),
  fgb = TRUE
)

# 渲染互動地圖至 RStudio Viewer
map_check <- mapview(
  comsa_ea_vis,
  zcol = "sero_status",
  col.regions = c("Selected (Serosurvey)" = "#FF3333", "Not Selected" = "#33CCFF"), # 鮮紅 vs 亮青藍
  color = "white",
  lwd = 1.0,
  alpha.regions = 0.35,
  label = comsa_ea_vis$tooltip_ea,
  layer.name = "COMSA EAs (341 Polygons)",
  popup = popupTable(
    comsa_ea_vis,
    zcol = c("cluster", "distrito", "posto", "localidade", "sero", "sero_status")
  )
) +
  mapview(
    summarized_moz_age_sf,
    col.regions = "#FFFF00", # 螢光亮黃實心點
    color = "black",         # 黑色輪廓增強立體感
    lwd = 1.2,
    cex = 6,                 # 醒目半徑
    alpha.regions = 0.95,
    label = summarized_moz_age_sf$tooltip_moz,
    layer.name = "Moz 30 Centroids (summarized_moz_age_sf)",
    popup = popupTable(
      summarized_moz_age_sf,
      zcol = c("csa_code", "n_households", "n_individuals", "seroprev", "mean_lon", "mean_lat")
    )
  )

print(map_check)
cat("  - 互動地圖已成功輸出至 Viewer 面板。\n")


# ------------------------------------------------------------------------------
# 步驟 4：成對大圓距離配對與 seroprev 轉錄 (Moz -> COMSA)
# ------------------------------------------------------------------------------
cat("\n>>> 執行步驟 4：透過最近鄰大圓距離配對轉錄 seroprev...\n")

# 切換 s2 開啟狀態以使用高精度球面大圓距離計算（免安裝 lwgeom）
sf_use_s2(TRUE)
dist_matrix <- st_distance(summarized_moz_age_sf, summarized_comsa_sf)
sf_use_s2(FALSE)

# 對每個 Moz cluster 找到距離最近的 COMSA cluster
nearest_comsa_indices <- apply(dist_matrix, 1, which.min)
nearest_distances_km  <- apply(dist_matrix, 1, min) / 1000

# 建立對照配對表
matching_df <- data.frame(
  csa_code = summarized_moz_age_sf$csa_code,
  matched_comsa_cluster = summarized_comsa_sf$cluster[nearest_comsa_indices],
  dist_km = nearest_distances_km,
  seroprev = summarized_moz_age_sf$seroprev
)

# 驗證配對品質
cat("  - 配對最大距離 (Max distance):", round(max(matching_df$dist_km), 2), "km\n")
cat("  - 配對中位數距離 (Median distance):", round(median(matching_df$dist_km), 2), "km\n")

if (length(unique(matching_df$matched_comsa_cluster)) != 30) {
  warning("警告：存在多個 Moz Cluster 配對到同一個 COMSA Cluster 的情況！")
} else {
  cat("  - 驗證成功：30 個 Moz Clusters 與 COMSA 達成嚴格 1 對 1 配對。\n")
}

# 將 seroprev 轉錄填入 COMSA 多邊形與 Centroid 物件
comsa_dissolved <- comsa_dissolved %>%
  left_join(
    matching_df %>% select(matched_comsa_cluster, seroprev),
    by = c("cluster" = "matched_comsa_cluster")
  )

summarized_comsa_sf <- summarized_comsa_sf %>%
  left_join(
    matching_df %>% select(matched_comsa_cluster, seroprev),
    by = c("cluster" = "matched_comsa_cluster")
  )


# ------------------------------------------------------------------------------
# 步驟 5：Base AC 抽樣母體點位標註與最近中心點裁決
# ------------------------------------------------------------------------------
cat("\n>>> 執行步驟 5：Base AC 4,874 點位空間相交檢驗與代表點裁決...\n")

# 讀取 Base AC 4,874 點位
base_ac_sf <- st_read(
  here("01_data", "raw", "ZBAC_Points_Final", "04_Base_AC_Final.shp"),
  quiet = TRUE
) %>% st_transform(4326)

# 建立臨時列編號以利後續 left_join
base_ac_sf$temp_row_id <- seq_len(nrow(base_ac_sf))

# 空間包含檢驗：找出實體落在 115 個 COMSA 多邊形內的所有候選點
pts_in_poly <- st_join(
  base_ac_sf,
  comsa_dissolved %>% select(comsa_cluster = cluster, poly_sero = sero, poly_seroprev = seroprev),
  join = st_intersects,
  left = FALSE
)

cat("  - 實體落入 COMSA 多邊形內的 Base AC 候選點總數:", nrow(pts_in_poly), "\n")

# 多點裁決機制：若一個 Cluster 多邊形內有多個點，精準計算至幾何中心的距離
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

# 每個 cluster 僅挑選距離幾何中心最近的 1 個點作為該 cluster 代表點
primary_selected_pts <- pts_with_dist %>%
  group_by(comsa_cluster) %>%
  slice_min(order_by = dist_to_centroid, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("  - 經過最近中心點裁決後，確認指派的 COMSA 代表點數量:", nrow(primary_selected_pts), "\n")

# 將裁決代表點整併回原始的 base_ac_sf (4,874 點)
base_ac_processed <- base_ac_sf %>%
  left_join(
    st_drop_geometry(primary_selected_pts) %>% 
      select(temp_row_id, comsa_cluster, poly_sero, poly_seroprev),
    by = "temp_row_id"
  ) %>%
  mutate(
    # 建立分類變數：0 (母體底層), 1 (COMSA only), 2 (COMSA + sero)
    COMSA = case_when(
      !is.na(comsa_cluster) & poly_sero == 1 ~ 2,
      !is.na(comsa_cluster) & poly_sero == 0 ~ 1,
      TRUE                                   ~ 0
    ),
    comsa_cluster_id = comsa_cluster,
    seroprev = if_else(COMSA == 2, poly_seroprev, NA_real_)
  ) %>%
  select(-temp_row_id, -comsa_cluster, -poly_sero, -poly_seroprev)


# ------------------------------------------------------------------------------
# 步驟 6：資料驗證與持久化儲存
# ------------------------------------------------------------------------------
cat("\n>>> 執行步驟 6：資料完整性驗證與 RDS 檔案匯出...\n")

# 還原 s2 球面幾何引擎設定
sf_use_s2(TRUE)

# 驗證總筆數與欄位分類分佈
cat("總筆數確認 (應精確為 4874):", nrow(base_ac_processed), "\n")
cat("COMSA 欄位分類統計：\n")
print(table(COMSA = base_ac_processed$COMSA, useNA = "always"))

cat("\n各類別具有 seroprev 數值的點位數：\n")
print(table(COMSA = base_ac_processed$COMSA, Has_Seroprev = !is.na(base_ac_processed$seroprev)))

# 確保輸出目錄存在並寫入 RDS
output_dir <- here("01_data", "processed")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_path <- here("01_data", "processed", "base_ac_processed.rds")
saveRDS(base_ac_processed, file = output_path)

cat("\n處理完成！檔案已安全封裝儲存至：\n", output_path, "\n")