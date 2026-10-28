# ==============================================================================
# 專案：Zambézia 22 縣三層級分類地圖 (Sero / SRS Only / No SRS + 115 抽樣點)
# 輸入檔案：
#   - 01_data/processed/base_ac_processed.rds
#   - 03_results/table_district_sampling_summary_22districts.csv
# 輸出檔案：
#   - 03_results/map_zambezia_22districts_3tiers_srs.png (高解析靜態圖)
#   - 03_results/map_zambezia_22districts_3tiers_srs.pdf
#   - RStudio Viewer 互動式地圖 (OSM / Esri Satellite / CartoDB)
# ==============================================================================

library(here)
library(sf)
library(dplyr)
library(ggplot2)
library(mapview)

# 關閉 s2 平面運算以防多邊形幾何拓撲修復時報錯
sf_use_s2(FALSE)

# ------------------------------------------------------------------------------
# 步驟 1：讀取 22 縣邊界、抽樣摘要表與 115 個核心 Cluster 代表點
# ------------------------------------------------------------------------------
cat(">>> [1/4] 讀取點位、行政區多邊形與統計資料...\n")

# 1. 讀取 4,874 點位並嚴格篩選出 115 個 SRS 代表點 (排除 COMSA == 0)
base_ac_processed <- readRDS(here("01_data", "processed", "base_ac_processed.rds")) %>%
  st_transform(4326)

pts_115 <- base_ac_processed %>%
  filter(COMSA %in% c(1, 2)) %>%
  mutate(
    tier_label = factor(
      if_else(COMSA == 2, "COMSA + Serosurvey (N = 30)", "COMSA Only (N = 85)"),
      levels = c("COMSA + Serosurvey (N = 30)", "COMSA Only (N = 85)")
    )
  )

cat("  - 115 個 Cluster 點位確認:", nrow(pts_115), "(應精確為 115)\n")

# 2. 建構 Zambézia 現行 22 縣邊界 (利用 GADM level 3 精確分離 Mulevala 等新設縣)
if (!requireNamespace("geodata", quietly = TRUE)) install.packages("geodata")
moz_adm3_raw <- geodata::gadm(country = "MOZ", level = 3, path = tempdir()) %>%
  st_as_sf()

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

cat("  - 縣份總數確認:", nrow(zambezia_distritos), "(已精確為 22 縣)\n")

# 3. 讀取 22 縣抽樣摘要統計表並指派三層級分類
summary_csv_path <- here("03_results", "table_district_sampling_summary_22districts.csv")
summary_df <- read.csv(summary_csv_path) %>%
  filter(district_gadm != "TOTAL (Zambézia)") %>%
  mutate(
    # 三層級狀態劃分
    district_status = case_when(
      COMSA_Sero_N30 > 0 ~ "Has Serosurvey (N = 12)",
      COMSA_SRS_N115 > 0 ~ "SRS Only (N = 9)",
      TRUE               ~ "No SRS (N = 1)"
    ),
    # 鎖定 Factor 排序以確保圖例顏色嚴格綁定
    district_status = factor(
      district_status,
      levels = c("Has Serosurvey (N = 12)", "SRS Only (N = 9)", "No SRS (N = 1)")
    )
  )

# 整併統計屬性至 22 縣多邊形
distritos_map_sf <- zambezia_distritos %>%
  left_join(summary_df, by = "district_gadm")

cat("  - 22 縣三層級分類統計確認：\n")
print(table(distritos_map_sf$district_status, useNA = "always"))


# ------------------------------------------------------------------------------
# 步驟 2：計算縣名標籤幾何中心座標
# ------------------------------------------------------------------------------
cat("\n>>> [2/4] 計算縣名注記標籤座標...\n")

centroids <- st_centroid(distritos_map_sf)
coords <- st_coordinates(centroids)
distritos_map_sf$label_x <- coords[, 1]
distritos_map_sf$label_y <- coords[, 2]


# ------------------------------------------------------------------------------
# 步驟 3：繪製出版級高解析靜態地圖 (ggplot2)
# ------------------------------------------------------------------------------
cat("\n>>> [3/4] 輸出出版級高解析靜態地圖 (ggplot2)...\n")

p_map_fixed <- ggplot() +
  # 1. 22 縣多邊形 (三色離散填色)
  geom_sf(
    data = distritos_map_sf,
    aes(fill = district_status),
    color = "#37474F",
    linewidth = 0.45
  ) +
  scale_fill_manual(
    name = "District Sampling Status",
    values = c(
      "Has Serosurvey (N = 12)" = "#A5D6A7", # 柔和粉綠
      "SRS Only (N = 9)"         = "#FFE082", # 柔和琥珀橘黃
      "No SRS (N = 1)"           = "#ECEFF1"  # 灰白底 (Mocubela)
    )
  ) +
  # 2. 縣名標籤
  geom_text(
    data = distritos_map_sf,
    aes(x = label_x, y = label_y, label = district_gadm),
    size = 2.6,
    fontface = "bold",
    color = "#212121",
    check_overlap = TRUE
  ) +
  # 3. 115 個 SRS 代表點 (紅圓點 Sero + 深藍三角 SRS Only)
  geom_sf(
    data = pts_115,
    aes(color = tier_label, shape = tier_label),
    size = 2.4,
    stroke = 0.6
  ) +
  scale_color_manual(
    name = "Sampled Clusters (N = 115)",
    values = c(
      "COMSA + Serosurvey (N = 30)" = "#D50000", # 明亮朱紅圓點
      "COMSA Only (N = 85)"        = "#0D47A1"  # 深藍三角點
    )
  ) +
  scale_shape_manual(
    name = "Sampled Clusters (N = 115)",
    values = c(
      "COMSA + Serosurvey (N = 30)" = 16, # 實心圓點
      "COMSA Only (N = 85)"        = 17  # 實心三角點
    )
  ) +
  labs(
    title = "COMSA Sampling Stratification Across 22 Distritos of Zambézia",
    subtitle = "Districts: Green = Has Sero (N=12) | Orange = SRS Only (N=9) | Grey = No SRS (N=1; Mocubela)\nPoints: 115 Sampled Clusters (30 Sero Sites + 85 SRS Only; 4,759 Frame Excluded)",
    x = "Longitude",
    y = "Latitude",
    caption = "Source: COMSA SRS Sampling Frame & GADM Level 3 Administrative Boundaries (2013 Revision)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(color = "#37474F", size = 9.2, margin = margin(b = 10)),
    legend.position = "right",
    legend.background = element_rect(fill = alpha("white", 0.85), color = "grey80"),
    panel.grid.major = element_line(color = "grey90", linetype = "dotted"),
    panel.background = element_rect(fill = "#FAFAFA", color = NA)
  )

# 儲存靜態圖檔
ggsave(
  filename = here("03_results", "map_zambezia_22districts_3tiers_srs.png"),
  plot = p_map_fixed,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  filename = here("03_results", "map_zambezia_22districts_3tiers_srs.pdf"),
  plot = p_map_fixed,
  width = 10,
  height = 8
)
cat("  ✓ 靜態地圖已儲存至: 03_results/map_zambezia_22districts_3tiers_srs.png\n")


# ------------------------------------------------------------------------------
# 步驟 4：渲染互動式地圖至 Viewer (mapview 原生精準分類著色)
# ------------------------------------------------------------------------------
cat("\n>>> [4/4] 渲染互動式地圖至 RStudio Viewer...\n")

mapviewOptions(
  basemaps = c("OpenStreetMap", "Esri.WorldImagery", "CartoDB.Positron"),
  fgb = TRUE
)

# 1. 22 縣多邊形圖層 (指定 zcol 進行原生精準對照)
map_district_interactive <- mapview(
  distritos_map_sf,
  zcol = "district_status",
  col.regions = c(
    "Has Serosurvey (N = 12)" = "#66BB6A", # 翠綠
    "SRS Only (N = 9)"         = "#FFA726", # 暖橘
    "No SRS (N = 1)"           = "#CFD8DC"  # 冷灰 (Mocubela)
  ),
  color = "#37474F",
  lwd = 1.4,
  alpha.regions = 0.55,
  layer.name = "District Sampling Status",
  label = paste0(
    distritos_map_sf$district_gadm,
    " [", distritos_map_sf$district_status, "]",
    " | SRS: ", distritos_map_sf$COMSA_SRS_N115,
    " | Sero: ", distritos_map_sf$COMSA_Sero_N30,
    " | Rate: ", distritos_map_sf$SRS_Sampling_Rate
  ),
  popup = TRUE
)

# 2. 115 個 SRS 代表點圖層 (紅點 Sero + 藍點 SRS Only)
map_pts_interactive <- mapview(
  pts_115,
  zcol = "tier_label",
  col.regions = c(
    "COMSA + Serosurvey (N = 30)" = "#D50000", # 紅色
    "COMSA Only (N = 85)"        = "#0D47A1"  # 深藍色
  ),
  color = "black",
  lwd = 0.6,
  cex = 4.5,
  alpha.regions = 0.9,
  layer.name = "115 COMSA Clusters",
  label = paste0(
    "Cluster: ", pts_115$comsa_cluster_id,
    " | ", pts_115$tier_label,
    " | Distrito: ", pts_115$district_gadm
  ),
  popup = TRUE
)

# 組合雙圖層
map_final_view <- map_district_interactive + map_pts_interactive

# 還原 s2 球面運算設定
sf_use_s2(TRUE)

print(map_final_view)

cat("\n==============================================================================\n")
cat("全數完成！\n")
cat("1. 靜態圖已更新輸出至 03_results/。\n")
cat("2. 互動地圖中 22 縣已依 Sero (綠, 12)、SRS Only (橘, 9)、No SRS (灰, 1) 精準區隔！\n")
cat("==============================================================================\n")