# ==============================================================================
# 專案：莫三比克 Zambézia 省空間抽樣整合、成對距離分析與出版級圖表匯出
# 產出標的：
#   - 01_data/processed/base_ac_processed.rds
#   - 03_results/pairwise_distance_density_vertical.png (.pdf)
#   - 03_results/pairwise_distance_boxplot_horizontal.png (.pdf)
# ==============================================================================

library(here)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(mapview)
library(leafpop)
if (!requireNamespace("ggbeeswarm", quietly = TRUE)) install.packages("ggbeeswarm")
library(ggbeeswarm)

# 確保必要目錄結構存在
dir.create(here("01_data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("03_results"), recursive = TRUE, showWarnings = FALSE)

# 全域幾何設定：關閉 s2 球面幾何引擎，改採 GEOS 平面幾何引擎以杜絕多邊形自交報錯
sf_use_s2(FALSE)

# ------------------------------------------------------------------------------
# 步驟 1：COMSA EA 多邊形拓撲修復、融合與幾何中心摘要 (N = 115)
# ------------------------------------------------------------------------------
cat(">>> [1/8] 執行 COMSA 多邊形修復、融合與 Centroid 提取...\n")

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

cat("    - COMSA 融合後 Cluster 總數:", nrow(comsa_dissolved), "\n")
cat("    - 幾何中心點數量:", nrow(summarized_comsa_sf), "\n")


# ------------------------------------------------------------------------------
# 步驟 2：Mozambique 血清資料家戶幾何與 seroprev 彙總 (N = 30)
# ------------------------------------------------------------------------------
cat("\n>>> [2/8] 計算 Moz Age 30 個群聚之幾何中心與 seroprev...\n")

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
if (is.na(target_pos_col)) target_pos_col <- "pos"

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

cat("    - Moz Age 彙總群聚數:", nrow(summarized_moz_age_sf), "\n")


# ------------------------------------------------------------------------------
# 步驟 3：互動式衛星地圖檢驗 (檢查 30 個 Moz Centroids 分佈)
# ------------------------------------------------------------------------------
cat("\n>>> [3/8] 產出互動式地圖至 Viewer 面板...\n")

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
  layer.name = "COMSA EAs (341 Polygons)",
  popup = popupTable(
    comsa_ea_vis,
    zcol = c("cluster", "distrito", "posto", "localidade", "sero", "sero_status")
  )
) +
  mapview(
    summarized_moz_age_sf,
    col.regions = "#FFFF00",
    color = "black",
    lwd = 1.2,
    cex = 6,
    alpha.regions = 0.95,
    label = summarized_moz_age_sf$tooltip_moz,
    layer.name = "Moz 30 Centroids",
    popup = popupTable(
      summarized_moz_age_sf,
      zcol = c("csa_code", "n_households", "n_individuals", "seroprev", "mean_lon", "mean_lat")
    )
  )

print(map_check)


# ------------------------------------------------------------------------------
# 步驟 4：大圓距離配對與 seroprev 轉錄 (Moz -> COMSA)
# ------------------------------------------------------------------------------
cat("\n>>> [4/8] 透過最近鄰大圓距離配對轉錄 seroprev...\n")

# 切換 s2 開啟狀態以使用高精度球面大圓距離計算
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

cat("    - 配對最大距離 (Max distance):", round(max(matching_df$dist_km), 2), "km\n")
cat("    - 配對中位數距離 (Median distance):", round(median(matching_df$dist_km), 2), "km\n")

if (length(unique(matching_df$matched_comsa_cluster)) != 30) {
  warning("警告：存在多個 Moz Cluster 配對到同一個 COMSA Cluster 的情況！")
} else {
  cat("    - 驗證成功：30 個 Moz Clusters 與 COMSA 達成嚴格 1 對 1 配對。\n")
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
# 步驟 5：Base AC 空間相交檢驗、最近中心點裁決與欄位建立
# ------------------------------------------------------------------------------
cat("\n>>> [5/8] Base AC 4,874 點位空間相交檢驗與代表點裁決...\n")

base_ac_sf <- st_read(
  here("01_data", "raw", "ZBAC_Points_Final", "04_Base_AC_Final.shp"),
  quiet = TRUE
) %>% st_transform(4326)

base_ac_sf$temp_row_id <- seq_len(nrow(base_ac_sf))

# 空間包含檢驗：找出實體落在 115 個 COMSA 多邊形內的所有候選點
pts_in_poly <- st_join(
  base_ac_sf,
  comsa_dissolved %>% select(comsa_cluster = cluster, poly_sero = sero, poly_seroprev = seroprev),
  join = st_intersects,
  left = FALSE
)

cat("    - 實體落入 COMSA 多邊形內的 Base AC 候選點總數:", nrow(pts_in_poly), "\n")

# 多點裁決機制：計算候選點到幾何中心距離，僅保留距離幾何中心最近的 1 個點
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

cat("    - 經過最近中心點裁決後，確認指派的 COMSA 代表點數量:", nrow(primary_selected_pts), "\n")

# 將裁決代表點整併回原始的 base_ac_sf (4,874 點)
base_ac_processed <- base_ac_sf %>%
  left_join(
    st_drop_geometry(primary_selected_pts) %>% 
      select(temp_row_id, comsa_cluster, poly_sero, poly_seroprev),
    by = "temp_row_id"
  ) %>%
  mutate(
    # 分類標籤：0 (母體底層), 1 (COMSA only 85), 2 (COMSA + sero 30)
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
# 步驟 6：Master 資料集儲存至 01_data/processed/
# ------------------------------------------------------------------------------
cat("\n>>> [6/8] 驗證資料完整性並寫入 RDS...\n")

output_rds_path <- here("01_data", "processed", "base_ac_processed.rds")
saveRDS(base_ac_processed, file = output_rds_path)

cat("    - 總筆數 (確認為 4,874):", nrow(base_ac_processed), "\n")
cat("    - COMSA 分類計數:\n")
print(table(COMSA = base_ac_processed$COMSA, useNA = "always"))
cat("    - Master 資料集已成功寫入:", output_rds_path, "\n")


# ------------------------------------------------------------------------------
# 步驟 7：計算三個抽樣層級的大圓成對距離 (Pairwise Distance in km)
# ------------------------------------------------------------------------------
cat("\n>>> [7/8] 計算三個階層的成對大圓距離矩陣...\n")

# 開啟 s2 確保大圓球面距離計算準確度
sf_use_s2(TRUE)

calc_pairwise_km <- function(sf_obj, group_label) {
  dist_mat <- st_distance(sf_obj)
  dists_km <- as.numeric(dist_mat[lower.tri(dist_mat)]) / 1000
  data.frame(
    distance_km = dists_km,
    dataset = group_label,
    stringsAsFactors = FALSE
  )
}

# 提取各層物件
sf_frame     <- base_ac_processed
sf_comsa_srs <- base_ac_processed %>% filter(COMSA %in% c(1, 2))
sf_comsa_sero<- base_ac_processed %>% filter(COMSA == 2)

dist_frame_raw <- calc_pairwise_km(sf_frame, "Sampling Frame (N = 4,874)")
dist_srs       <- calc_pairwise_km(sf_comsa_srs, "COMSA SRS (N = 115)")
dist_sero      <- calc_pairwise_km(sf_comsa_sero, "COMSA Sero (N = 30)")

# 針對 Sampling frame 距離抽樣 30 萬筆以平衡視覺化記憶體負擔
set.seed(42)
if (nrow(dist_frame_raw) > 300000) {
  dist_frame <- dist_frame_raw[sample(nrow(dist_frame_raw), 300000), ]
} else {
  dist_frame <- dist_frame_raw
}

levels_order <- c(
  "Sampling Frame (N = 4,874)",
  "COMSA SRS (N = 115)",
  "COMSA Sero (N = 30)"
)

pairwise_df <- bind_rows(dist_frame, dist_srs, dist_sero) %>%
  mutate(dataset = factor(dataset, levels = levels_order))

palette_colors <- c(
  "Sampling Frame (N = 4,874)" = "#2B8CBE", # 海洋藍
  "COMSA SRS (N = 115)"        = "#E34A33", # 朱紅
  "COMSA Sero (N = 30)"        = "#2CA25F"  # 翠綠
)


# ------------------------------------------------------------------------------
# 步驟 8：繪製兩張圖形並輸出至 03_results/ (PNG + PDF)
# ------------------------------------------------------------------------------
cat("\n>>> [8/8] 繪製成對距離圖表並輸出檔案至 03_results/...\n")

# --- 圖 1：由上至下垂直核密度圖 ---
p1_density_vertical <- ggplot(pairwise_df, aes(x = distance_km, fill = dataset, color = dataset)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  facet_wrap(~ dataset, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = palette_colors) +
  scale_color_manual(values = palette_colors) +
  scale_x_continuous(limits = c(0, max(pairwise_df$distance_km) * 1.05), expand = c(0, 0)) +
  labs(
    title = "Pairwise Geodesic Distance: Density Distributions",
    subtitle = "Vertical facet comparison across sampling hierarchy in Zambézia",
    x = "Pairwise Distance (km)",
    y = "Kernel Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray93", color = NA),
    plot.title = element_text(face = "bold", size = 14),
    panel.spacing = unit(1.0, "lines"),
    axis.title = element_text(face = "bold")
  )

# --- 圖 2：由左至右水平排列的資料點分布 + 箱形圖 ---
set.seed(42)
sampled_points_df <- pairwise_df %>%
  group_by(dataset) %>%
  slice_sample(n = 2500) %>%
  ungroup()

p2_box_beeswarm_horizontal <- ggplot(pairwise_df, aes(x = dataset, y = distance_km, fill = dataset)) +
  # --- 將原本的 geom_jitter 改為 geom_quasirandom ---
  geom_quasirandom(
    data = sampled_points_df,
    aes(color = dataset),
    width = 0.28,             # 蜂群橫向展開的最大寬度
    size = 0.7,               # 點的大小
    alpha = 0.30,             # 透明度
    varwidth = TRUE,          # 密度越高的地方蜂群展開越寬（呈現輪廓形狀）
    show.legend = FALSE
  ) +
  # 上層維持半透明箱形圖
  geom_boxplot(
    width = 0.35,
    alpha = 0.75,
    color = "gray20",
    outlier.shape = NA,
    linewidth = 0.7
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 3.2,
    fill = "white",
    color = "black",
    stroke = 1.0
  ) +
  scale_fill_manual(values = palette_colors) +
  scale_color_manual(values = palette_colors) +
  labs(
    title = "Pairwise Distance Distribution: Beeswarm & Box Plots",
    subtitle = "Horizontal comparison (Diamonds denote group means; points show beeswarm distribution)",
    x = "Sampling Tier",
    y = "Pairwise Distance (km)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

# --- 檔案匯出至 03_results/ (PNG 300 DPI + 向量 PDF) ---

# 圖 1 儲存
ggsave(
  filename = here("03_results", "pairwise_distance_density_vertical.png"),
  plot = p1_density_vertical,
  width = 8.5, height = 7.5, dpi = 300
)
ggsave(
  filename = here("03_results", "pairwise_distance_density_vertical.pdf"),
  plot = p1_density_vertical,
  width = 8.5, height = 7.5
)

# 圖 2 儲存
ggsave(
  filename = here("03_results", "pairwise_distance_boxplot_horizontal.png"),
  plot = p2_box_beeswarm_horizontal,
  width = 9.0, height = 6.0, dpi = 300
)
ggsave(
  filename = here("03_results", "pairwise_distance_boxplot_horizontal.pdf"),
  plot = p2_box_beeswarm_horizontal,
  width = 9.0, height = 6.0
)

cat("\n==============================================================================\n")
cat("全流程順利執行完畢！產出檢查：\n")
cat("1. RDS 檔案：", output_rds_path, "\n")
cat("2. 圖 1 (PNG/PDF)：03_results/pairwise_distance_density_vertical.[png|pdf]\n")
cat("3. 圖 2 (PNG/PDF)：03_results/pairwise_distance_boxplot_horizontal.[png|pdf]\n")
cat("==============================================================================\n")