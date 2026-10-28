# 載入必要套件
library(here)
library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(mapview)
library(leafpop)

# 1. 驗證專案根目錄基準點
print(here())

# 2. 讀取莫三比克 Zambézia 省邊界作為底圖
moz_adm1 <- st_read(here("01_data", "raw", "moz_adm_20190607b_shp", "moz_admbnda_adm1_ine_20190607.shp"))
zambezia_border <- moz_adm1 %>% 
  filter(ADM1_PT == "Zambezia") %>% 
  st_transform(4326)

# 3. 讀取 Shapefile 資料
# comsa_sf 讀入後為 MULTIPOLYGON，需計算 centroid 轉為代表點
comsa_raw   <- st_read(here("01_data", "raw", "comsa.zambezia_with_clusters", "comsa.zambezia_with_clusters.shp"))
comsa_sf    <- st_centroid(st_transform(comsa_raw, 4326))

# base_ac_sf 為點位圖層 (POINT)
base_ac_sf  <- st_read(here("01_data", "raw", "ZBAC_Points_Final", "04_Base_AC_Final.shp")) %>% 
  st_transform(4326)

# 4. 讀取 RDS 資料檔並轉換為 sf (使用全小寫欄位 lon, lat)
moz_age <- readRDS(here("01_data", "raw", "moz_age.rds"))
moz_age_sf <- st_as_sf(moz_age, coords = c("lon", "lat"), crs = 4326, remove = FALSE)


# 確保所有圖層 CRS 一致
common_crs <- st_crs(4326)
comsa_sf        <- st_transform(comsa_sf, common_crs)
base_ac_sf      <- st_transform(base_ac_sf, common_crs)
moz_age_sf      <- st_transform(moz_age_sf, common_crs)
zambezia_border <- st_transform(zambezia_border, common_crs)

# --- 步驟 5：繪製個別分佈圖 ---

p_base <- ggplot() +
  geom_sf(data = zambezia_border, fill = "gray96", color = "gray60", linewidth = 0.4) +
  theme_minimal(base_size = 11) +
  labs(x = "Longitude", y = "Latitude")

p1 <- p_base + 
  geom_sf(data = moz_age_sf, color = "#2CA25F", size = 0.8, alpha = 0.5) +
  labs(title = paste0("Moz Age (N = ", nrow(moz_age_sf), ")"))

p2 <- p_base + 
  geom_sf(data = comsa_sf, color = "#E34A33", size = 1.5, alpha = 0.7) +
  labs(title = paste0("COMSA Centroids (N = ", nrow(comsa_sf), ")"))

p3 <- p_base + 
  geom_sf(data = base_ac_sf, color = "#2B8CBE", size = 1.2, alpha = 0.6) +
  labs(title = paste0("Base AC Final (N = ", nrow(base_ac_sf), ")"))

# 橫向並排三張子圖
p_separate <- p1 | p2 | p3
print(p_separate)

# --- 步驟 6：繪製重疊圖（檢查點位是否重合） ---

p_overlay <- p_base +
  geom_sf(data = moz_age_sf, aes(color = "Moz Age"), size = 0.9, alpha = 0.4) +
  geom_sf(data = comsa_sf, aes(color = "COMSA Centroids"), size = 1.8, shape = 17, alpha = 0.8) +
  geom_sf(data = base_ac_sf, aes(color = "Base AC"), size = 1.5, shape = 15, alpha = 0.8) +
  scale_color_manual(
    name = "Dataset",
    values = c("Moz Age" = "#2CA25F", "COMSA Centroids" = "#E34A33", "Base AC" = "#2B8CBE")
  ) +
  labs(
    title = "Spatial Overlay of Sampling Datasets in Zambézia",
    subtitle = "Checking spatial alignment across datasets"
  ) +
  theme(legend.position = "right")

print(p_overlay)

# --- 步驟 7：互動式地圖（縮放精確確認重合度） ---

mapview(zambezia_border, alpha.regions = 0.05, color = "black", layer.name = "Zambézia Border") +
  mapview(moz_age_sf, col.regions = "#2CA25F", cex = 3, layer.name = "Moz Age") +
  mapview(comsa_sf, col.regions = "#E34A33", cex = 5, layer.name = "COMSA Centroids") +
  mapview(base_ac_sf, col.regions = "#2B8CBE", cex = 4, layer.name = "Base AC")

# 檢查各圖層點數
cat("--- 資料集點位數量統計 ---\n")
cat("Moz age 點數 (N):", nrow(moz_age_sf), "\n")
cat("COMSA centroids 點數 (N):", nrow(comsa_sf), "\n")
cat("Base AC 點數 (N):", nrow(base_ac_sf), "\n")



# 檢查不重複的 GPS 點數 (Unique coordinates)
cat("Moz age 不重複座標數:", nrow(distinct(st_drop_geometry(moz_age_sf), lon, lat)), "\n")
cat("COMSA 不重複群聚數:", nrow(distinct(st_drop_geometry(comsa_sf), cluster)), "\n")
cat("Base AC 不重複座標數:", length(unique(st_geometry(base_ac_sf))), "\n")


### Check EA for sero or not
# ------------------------------------------------------------------------------
# 1. 準備 341 個 EA 資料，建立視覺化分類標籤
# ------------------------------------------------------------------------------

# 讀取並轉換座標系統至 WGS84 (EPSG:4326)
comsa_ea_341 <- comsa_raw %>%
  st_transform(4326) %>%
  mutate(
    # 將數值型 sero 轉為清晰易讀的分類因子
    sero_status = case_when(
      sero == 1 ~ "Selected (Serosurvey)",
      sero == 0 ~ "Not Selected",
      TRUE      ~ "Unknown"
    ),
    sero_status = factor(sero_status, levels = c("Selected (Serosurvey)", "Not Selected"))
  )

# 統計被選中與未選中的 EA 數量
cat("--- 341 個 EA 的 Serosurvey 選取統計 ---\n")
print(table(comsa_ea_341$sero_status))

# ------------------------------------------------------------------------------
# 2. 靜態全省分佈圖 (ggplot2)
# ------------------------------------------------------------------------------

p_ea_sero <- ggplot() +
  # Zambézia 省界底圖
  geom_sf(data = zambezia_border, fill = "gray96", color = "gray60", linewidth = 0.5) +
  # 341 個 EA 多邊形：依 sero_status 著色
  geom_sf(
    data = comsa_ea_341,
    aes(fill = sero_status, color = sero_status),
    alpha = 0.65,
    linewidth = 0.3
  ) +
  # 設定對比配色：被選中設為深珊瑚紅 (#D95F02)，未選中設為冷灰藍 (#7570B3)
  scale_fill_manual(
    name = "Serosurvey Status",
    values = c("Selected (Serosurvey)" = "#D95F02", "Not Selected" = "#7570B3")
  ) +
  scale_color_manual(
    name = "Serosurvey Status",
    values = c("Selected (Serosurvey)" = "#A63603", "Not Selected" = "#542788")
  ) +
  labs(
    title = "COMSA Enumeration Areas (EAs) in Zambézia",
    subtitle = paste0("Total EAs = ", nrow(comsa_ea_341), " | Stratified by Serosurvey Selection (sero)"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "gray85"),
    panel.grid.minor = element_blank()
  )

print(p_ea_sero)

# ------------------------------------------------------------------------------
# 3. 互動式地圖 (mapview，免 API Key 圖磚)
# ------------------------------------------------------------------------------

# 設定使用免 API Key 的免費圖磚來源
mapviewOptions(basemaps = c("OpenStreetMap", "Esri.WorldImagery"))

mapview(
  zambezia_border,
  alpha.regions = 0.05,
  color = "black",
  lwd = 1.5,
  layer.name = "省界 Border"
) +
  mapview(
    comsa_ea_341,
    zcol = "sero_status",
    col.regions = c("#D95F02", "#7570B3"),
    alpha.regions = 0.6,
    lwd = 0.8,
    layer.name = "COMSA EAs (341 Polygons)",
    popup = popupTable(
      comsa_ea_341,
      zcol = c("cluster", "distrito", "posto", "localidade", "sero", "sero_status")
    )
  )



### Check where are 999 Polygons?

# 篩選出 999 的多邊形並統一投影至 WGS84
comsa_999 <- comsa_raw %>%
  filter(cluster == 999) %>%
  st_transform(4326)

# 設定免費開源圖磚
mapviewOptions(basemaps = c("OpenStreetMap", "Esri.WorldImagery"))

# 繪製地圖：以亮紅色突顯 999 的 5 個多邊形
mapview(
  comsa_999,
  col.regions = "red",
  color = "darkred",
  alpha.regions = 0.5,
  lwd = 2,
  layer.name = "Cluster 999 Polygons"
)


############# Inconsistent Sero for EAs within the same cluster #################

inconsistent_clusters <- comsa_raw %>%
  st_drop_geometry() %>%
  group_by(cluster) %>%
  summarize(
    n_ea = n(),
    n_distinct_sero = n_distinct(sero),
    sero_vals = paste(sort(unique(sero)), collapse = ", "),
    .groups = "drop"
  )

inconsistent_clusters_filtered <- inconsistent_clusters %>%
  filter(n_distinct_sero>1)

# 印出檢查結果
if (nrow(inconsistent_clusters_filtered) == 0) {
  cat("確認完畢：所有 cluster 內部的 EA 之 sero 標籤完全一致（沒有任何 cluster 同時存在 1 與 0）。\n")
} else {
  cat("注意：發現有", nrow(inconsistent_clusters_filtered), "個 cluster 內部的 EA sero 標籤不一致：\n")
  print(inconsistent_clusters_filtered)
}

sum(as.numeric(inconsistent_clusters$sero_vals))


########## Check If Neiboring EAs always belong to the same Cluster ##############


library(sf)
library(dplyr)
library(mapview)
library(leafpop)

if (!requireNamespace("spdep", quietly = TRUE)) install.packages("spdep")
library(spdep)

# ------------------------------------------------------------------------------
# 1. 聚合為 Cluster 多邊形並修復拓撲
# ------------------------------------------------------------------------------

sf_use_s2(FALSE) # 避免球面邊界自交問題

# 溶解同一個 cluster 的碎片，產出 116 個獨立 Cluster 多邊形
comsa_cluster_poly <- comsa_raw %>%
  st_transform(4326) %>%
  st_make_valid() %>%
  group_by(cluster) %>%
  summarize(
    sero = max(sero, na.rm = TRUE),
    geometry = st_union(st_combine(geometry)),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 2. 建立 Cluster 相鄰關係並進行圖著色 (Graph Coloring)
# ------------------------------------------------------------------------------

# 建立 116 個 Cluster 之間的空間相鄰關係 (Queen adjacency)
nb_cluster <- poly2nb(st_as_sf(comsa_cluster_poly))

# 定義 8 種肉眼辨識度極高且彼此反差強烈的顏色池 (Qualitative Palette)
cluster_color_pool <- c(
  "#E41A1C", # 鮮紅
  "#377EB8", # 寶藍
  "#4DAF4A", # 亮綠
  "#984EA3", # 深紫
  "#FF7F00", # 暖橘
  "#FFFF33", # 明黃
  "#A65628", # 赭褐
  "#F781BF"  # 亮粉
)

# 貪婪圖著色演算法：確保相鄰的 Cluster 絕不分到同色
n_clust <- nrow(comsa_cluster_poly)
assigned_color_id <- integer(n_clust)

for (i in seq_len(n_clust)) {
  # 找出第 i 個 cluster 的所有鄰居
  neighbors <- nb_cluster[[i]]
  # 取出鄰居已經用過的顏色
  used_colors <- assigned_color_id[neighbors[neighbors > 0 & neighbors < i]]
  # 從顏色池中選擇第一個沒被鄰居使用的顏色
  available_color <- which(!seq_along(cluster_color_pool) %in% used_colors)[1]
  assigned_color_id[i] <- available_color
}

# 將著色結果寫入 Cluster 表格
comsa_cluster_poly$cluster_color <- cluster_color_pool[assigned_color_id]

# ------------------------------------------------------------------------------
# 3. 將著色結果 Mapping 回原始的 341 個 EA 多邊形
# ------------------------------------------------------------------------------

comsa_ea_mapped <- comsa_raw %>%
  st_transform(4326) %>%
  left_join(
    st_drop_geometry(comsa_cluster_poly) %>% select(cluster, cluster_color),
    by = "cluster"
  ) %>%
  mutate(
    tooltip_label = paste0("Cluster: ", cluster, " (EA row: ", row, ", sero: ", sero, ")")
  )

sf_use_s2(TRUE) # 還原球面引擎

# ------------------------------------------------------------------------------
# 4. 繪製互動地圖 (mapview)
# ------------------------------------------------------------------------------

mapviewOptions(basemaps = c("OpenStreetMap", "Esri.WorldImagery"))

mapview(
  comsa_ea_mapped,
  col.regions = comsa_ea_mapped$cluster_color, # 依 cluster 著色
  color = "black",                            # 黑色外框呈現各 EA 切割邊界
  lwd = 1.2,
  alpha.regions = 0.65,
  legend = FALSE,
  label = comsa_ea_mapped$tooltip_label,      # 滑鼠懸停顯示 Cluster 與 EA 編號
  layer.name = "COMSA EAs (Cluster Topology Colored)",
  popup = popupTable(
    comsa_ea_mapped,
    zcol = c("cluster", "distrito", "posto", "localidade", "sero")
  )
)


#################################################################################
#-------------------- Overlay base_ac_sf and comsa_sf --------------------------#
#################################################################################
library(sf)
library(dplyr)
library(mapview)
library(leafpop)

# ------------------------------------------------------------------------------
# 1. 準備圖層資料並確保 CRS 統一為 EPSG:4326
# ------------------------------------------------------------------------------

# A. 省界底圖 (排除非法日期欄位，消除 GDAL Error 1 警告)
zambezia_border_clean <- zambezia_border %>% 
  st_transform(4326) %>% 
  select(-any_of(c("date", "validOn", "validTo")))

# B. COMSA 341 個 EA 多邊形：依 sero 建立清晰的視覺分類標籤
comsa_ea_sero <- comsa_raw %>%
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

# C. Base AC 4,874 個空間取樣點
base_ac_clean <- base_ac_sf %>%
  st_transform(4326)

# ------------------------------------------------------------------------------
# 2. 繪製互動式疊加地圖 (mapview)
# ------------------------------------------------------------------------------

# 設定免 API Key 的免費圖磚來源
mapviewOptions(
  basemaps = c("OpenStreetMap", "Esri.WorldImagery"),
  fgb = TRUE  # 啟用 FlatGeobuf 二進制傳輸加速 4800+ 點位渲染
)

# 疊圖順序：省界 -> EA 多邊形 -> 4874 個點位
mapview(
  zambezia_border_clean,
  alpha.regions = 0.03,
  color = "gray30",
  lwd = 1.2,
  layer.name = "省界 Border"
) +
  mapview(
    comsa_ea_sero,
    zcol = "sero_status",
    col.regions = c("Selected (Serosurvey)" = "#E6550D", "Not Selected" = "#756BB1"), # 橘紅 vs 灰紫
    color = "black",                                                                    # 黑色 EA 細邊界
    lwd = 0.8,
    alpha.regions = 0.45,
    label = comsa_ea_sero$tooltip_ea,
    layer.name = "COMSA EAs (341 Polygons)",
    popup = popupTable(
      comsa_ea_sero,
      zcol = c("cluster", "distrito", "posto", "localidade", "sero", "sero_status")
    )
  ) +
  mapview(
    base_ac_clean,
    col.regions = "#004D80", # 深青藍實心點
    color = "white",         # 白色細邊圈，在深色與淺色背景上皆具高辨識度
    lwd = 0.5,
    cex = 3,                 # 點位半徑大小（適中不遮蔽多邊形）
    alpha.regions = 0.75,
    layer.name = "Base AC Points (N=4874)"
  )

#################################################################################
#-------------------- Overlay moz_age_sf and comsa_sf --------------------------#
#################################################################################


library(sf)
library(dplyr)
library(mapview)
library(leafpop)

# ------------------------------------------------------------------------------
# 1. 準備圖層資料並確保 CRS 統一為 EPSG:4326
# ------------------------------------------------------------------------------

# A. 省界底圖 (排除非法日期欄位，消除 GDAL Error 1 警告)
zambezia_border_clean <- zambezia_border %>% 
  st_transform(4326) %>% 
  select(-any_of(c("date", "validOn", "validTo")))

# B. COMSA 341 個 EA 多邊形：依 sero 區分顏色
comsa_ea_sero <- comsa_raw %>%
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

# C. moz_age 住戶點位：依 hhid 與不重複座標聚合 (N = 1,136 households)
moz_hh_sf <- moz_age %>%
  group_by(hhid, cluster, csa_code, district) %>%
  summarize(
    n_individuals = n(),               # 該住戶參與檢測的人數
    lon = first(lon),
    lat = first(lat),
    .groups = "drop"
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  mutate(
    tooltip_hh = paste0("HHID: ", hhid, " | Cluster: ", cluster, " (", n_individuals, " people)")
  )

cat("--- 圖層驗證 ---\n")
cat("COMSA EAs (N):", nrow(comsa_ea_sero), "\n")
cat("Moz Age Households (N):", nrow(moz_hh_sf), "\n")

# ------------------------------------------------------------------------------
# 2. 繪製互動式地圖 (mapview)
# ------------------------------------------------------------------------------

# 設定免 API Key 的免費圖磚
mapviewOptions(
  basemaps = c("OpenStreetMap", "Esri.WorldImagery"),
  fgb = TRUE
)

# 疊圖順序：省界 -> EA 多邊形 (sero 著色) -> Moz Age 家戶點
mapview(
  zambezia_border_clean,
  alpha.regions = 0.03,
  color = "gray30",
  lwd = 1.2,
  layer.name = "省界 Border"
) +
  mapview(
    comsa_ea_sero,
    zcol = "sero_status",
    col.regions = c("Selected (Serosurvey)" = "#E6550D", "Not Selected" = "#756BB1"), # 橘紅 vs 灰紫
    color = "black",
    lwd = 0.8,
    alpha.regions = 0.45,
    label = comsa_ea_sero$tooltip_ea,
    layer.name = "COMSA EAs (341 Polygons)",
    popup = popupTable(
      comsa_ea_sero,
      zcol = c("cluster", "distrito", "posto", "localidade", "sero", "sero_status")
    )
  ) +
  mapview(
    moz_hh_sf,
    col.regions = "#00E676",  # 高對比亮綠色實心點
    color = "#004D20",        # 深綠色外框線，強化在亮色底圖上的辨識度
    lwd = 0.8,
    cex = 4,                  # 點位半徑適中
    alpha.regions = 0.85,
    label = moz_hh_sf$tooltip_hh,
    layer.name = "Moz Age Households (N=1136)",
    popup = popupTable(
      moz_hh_sf,
      zcol = c("hhid", "cluster", "csa_code", "district", "n_individuals", "lon", "lat")
    )
  )


##### Check Topology

library(sf)
library(dplyr)
library(mapview)
library(leafpop)

# ------------------------------------------------------------------------------
# 1. 準備資料圖層 (EPSG:4326)
# ------------------------------------------------------------------------------

# A. 省界底圖 (排除非法日期欄位)
zambezia_border_clean <- zambezia_border %>% 
  st_transform(4326) %>% 
  select(-any_of(c("date", "validOn", "validTo")))

# B. COMSA 341 個 EA 多邊形：改用亮螢光色框線，避免遮蔽底層房屋與樹木
comsa_ea_sero <- comsa_raw %>%
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

# C. Moz Age 住戶點 (1,136 個家戶)
moz_hh_sf <- moz_age %>%
  group_by(hhid, cluster, csa_code, district) %>%
  summarize(
    n_individuals = n(),
    lon = first(lon),
    lat = first(lat),
    .groups = "drop"
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  mutate(
    tooltip_hh = paste0("HHID: ", hhid, " | Cluster: ", cluster, " (", n_individuals, " people)")
  )

# ------------------------------------------------------------------------------
# 2. 設定高解析度衛星與立體地形底圖
# ------------------------------------------------------------------------------

# 將 Esri.WorldImagery 設在第 1 位（預設載入即為高解析衛星空照圖）
mapviewOptions(
  basemaps = c(
    "Esri.WorldImagery",       # 高解析真實衛星影像（可辨識房屋結構、樹冠、農田）
    "Esri.WorldShadedRelief",  # 立體山脈陰影地形圖（看清山脈、山谷與丘陵起伏）
    "OpenStreetMap"            # 道路聚落向量底圖
  ),
  fgb = TRUE
)

# ------------------------------------------------------------------------------
# 3. 繪製互動地圖 (透明化填色以完整透出衛星影像細節)
# ------------------------------------------------------------------------------

mapview(
  zambezia_border_clean,
  alpha.regions = 0,         # 內部全透明，不遮擋衛星影像
  color = "yellow",          # 醒目亮黃色省界線
  lwd = 1.8,
  layer.name = "省界 Border"
) +
  mapview(
    comsa_ea_sero,
    zcol = "sero_status",
    col.regions = c("Selected (Serosurvey)" = "#FF3333", "Not Selected" = "#33CCFF"), # 鮮紅 vs 亮青天藍
    color = "white",           # 白色多邊形外框線
    lwd = 1.0,
    alpha.regions = 0.20,      # 多邊形降至 20% 透明度，使地面的房屋與樹木清晰可見
    label = comsa_ea_sero$tooltip_ea,
    layer.name = "COMSA EAs (341 Polygons)",
    popup = popupTable(
      comsa_ea_sero,
      zcol = c("cluster", "distrito", "posto", "localidade", "sero", "sero_status")
    )
  ) +
  mapview(
    moz_hh_sf,
    col.regions = "#FFFF00",   # 螢光亮黃實心點（在深綠色樹林與屋頂上對比最強）
    color = "black",           # 黑色細外框增加輪廓立體感
    lwd = 1.0,
    cex = 4.5,                 # 點位大小適中
    alpha.regions = 0.9,
    label = moz_hh_sf$tooltip_hh,
    layer.name = "Moz Age Households (N=1136)",
    popup = popupTable(
      moz_hh_sf,
      zcol = c("hhid", "cluster", "csa_code", "district", "n_individuals", "lon", "lat")
    )
  )






library(sf)
library(dplyr)
library(ggplot2)
library(mapview)

# --- 1. 確保多邊形與聚合中心點資料準備就緒 ---

# 原始多邊形圖層 (排除 999 雜訊，共 336 個多邊形)
comsa_poly_clean <- comsa_raw %>%
  st_transform(4326) %>%
  filter(cluster != 999)

# 聚合後的 115 個群聚中心點 (避開幾何自交報錯)
sf_use_s2(FALSE)  # 關閉 s2 避免多邊形自交報錯
comsa_sf_clean <- comsa_poly_clean %>%
  st_make_valid() %>%
  group_by(cluster) %>%
  summarize(.groups = "drop") %>%
  st_centroid()
sf_use_s2(TRUE)   # 重新開啟球面引擎

# --- 2. 靜態全區疊合圖 (ggplot2) ---

p_check_centroids <- ggplot() +
  # Zambézia 省界底圖
  geom_sf(data = zambezia_border, fill = "gray97", color = "gray60", linewidth = 0.4) +
  # 原始多邊形 (淡紅色半透明填充 + 細邊框)
  geom_sf(data = comsa_poly_clean, fill = "#FCAE91", color = "#FB6A4A", alpha = 0.45, linewidth = 0.25) +
  # 聚合中心點 (深藍色實心點 + 白色外圈增加對比)
  geom_sf(data = comsa_sf_clean, color = "white", size = 2.4, shape = 19) +
  geom_sf(data = comsa_sf_clean, color = "#08519C", size = 1.6, shape = 19) +
  labs(
    title = "COMSA Raw Polygons vs. Aggregated Cluster Centroids",
    subtitle = paste0("Polygons (N = ", nrow(comsa_poly_clean), ") overlaid with Centroids (N = ", nrow(comsa_sf_clean), ")"),
    caption = "Light Red: Raw Polygon Patches | Blue Dots: Summarized Cluster Centroids",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "gray30")
  )

print(p_check_centroids)

# --- 3. 互動式地圖檢驗 (推薦：可放大至微觀層級檢查) ---

# 設定互動圖層
mapview(zambezia_border, 
        alpha.regions = 0.05, 
        color = "black", 
        lwd = 1.5, 
        layer.name = "省界 Border") +
  mapview(comsa_poly_clean, 
          zcol = "cluster", 
          legend = FALSE, 
          alpha.regions = 0.35, 
          layer.name = "原始多邊形 (Raw Polygons)") +
  mapview(comsa_sf_clean, 
          col.regions = "blue", 
          color = "white", 
          cex = 5, 
          layer.name = "聚合中心點 (Centroids)")


# ==============================================================================
# 取得 moz_age 30 個 Cluster 的中心點 (lat & lon) 並繪圖
# ==============================================================================

# 1. 計算 30 個 cluster 的平均經緯度 (依 csa_code 或 cluster 分組)
moz_clusters_30 <- moz_age %>%
  group_by(cluster, csa_code, district) %>%
  summarize(
    n_samples = n(),                 # 該群聚的樣本數
    lon = mean(lon, na.rm = TRUE),   # 群聚中心經度
    lat = mean(lat, na.rm = TRUE),   # 群聚中心緯度
    .groups = "drop"
  )

# 檢查確認剛好 30 筆
cat("Moz 30 Clusters 筆數:", nrow(moz_clusters_30), "\n")
print(head(moz_clusters_30))

# 2. 轉換為 sf 空間點位物件 (EPSG:4326)
moz_clusters_30_sf <- st_as_sf(moz_clusters_30, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# ------------------------------------------------------------------------------
# 3. 繪圖：三圖獨立呈現 (30 Clusters vs 116 COMSA vs 4874 Base AC)
# ------------------------------------------------------------------------------
p_base <- ggplot() +
  geom_sf(data = zambezia_border, fill = "gray96", color = "gray60", linewidth = 0.4) +
  theme_minimal(base_size = 11) +
  labs(x = "Longitude", y = "Latitude")

# A. Moz 30 Clusters 圖
p_moz30 <- p_base +
  geom_sf(data = moz_clusters_30_sf, color = "#2CA25F", size = 2.2, alpha = 0.8) +
  labs(title = paste0("Moz Age Clusters (N = ", nrow(moz_clusters_30_sf), ")"))

# B. COMSA Centroids (N = 116)
p_comsa <- p_base +
  geom_sf(data = comsa_sf, color = "#E34A33", size = 1.6, alpha = 0.7) +
  labs(title = paste0("COMSA Clusters (N = ", nrow(comsa_sf), ")"))

# C. Base AC Points (N = 4874)
p_baseac <- p_base +
  geom_sf(data = base_ac_sf, color = "#2B8CBE", size = 0.8, alpha = 0.5) +
  labs(title = paste0("Base AC Points (N = ", nrow(base_ac_sf), ")"))

# 橫向並排三張子圖
p_clusters_separate <- p_moz30 | p_comsa | p_baseac
print(p_clusters_separate)

# ------------------------------------------------------------------------------
# 4. 繪圖：疊加圖（檢查 30 個群聚是否剛好落在 COMSA 與 Base AC 上）
# ------------------------------------------------------------------------------
p_overlay_30 <- p_base +
  geom_sf(data = base_ac_sf, aes(color = "Base AC (N=4874)"), size = 0.8, alpha = 0.3) +
  geom_sf(data = comsa_sf, aes(color = "COMSA Clusters (N=116)"), size = 2, shape = 17, alpha = 0.7) +
  geom_sf(data = moz_clusters_30_sf, aes(color = "Moz Clusters (N=30)"), size = 3, shape = 19, stroke = 0.8) +
  scale_color_manual(
    name = "Dataset",
    values = c(
      "Moz Clusters (N=30)"     = "#2CA25F", 
      "COMSA Clusters (N=116)"  = "#E34A33", 
      "Base AC (N=4874)"        = "#2B8CBE"
    )
  ) +
  labs(
    title = "Spatial Overlay: 30 Moz Clusters vs COMSA vs Base AC",
    subtitle = "Zambézia Province, Mozambique"
  ) +
  theme(legend.position = "right")

print(p_overlay_30)

# ------------------------------------------------------------------------------
# 5. 互動式地圖（用 mapview 放大確認 30 個點與其他點位是否完美重合）
# ------------------------------------------------------------------------------
mapview(zambezia_border, alpha.regions = 0.05, color = "black", layer.name = "Zambézia Border") +
  mapview(moz_clusters_30_sf, col.regions = "#2CA25F", cex = 6, layer.name = "Moz Clusters (N=30)") +
  mapview(comsa_sf, col.regions = "#E34A33", cex = 4, layer.name = "COMSA (N=116)") +
  mapview(base_ac_sf, col.regions = "#2B8CBE", cex = 2, layer.name = "Base AC (N=4874)")


# ==============================================================================
# 取得pairwise distance並繪圖
# ==============================================================================

library(sf)
library(dplyr)
library(ggplot2)

# --- 1. 定義計算成對大圓距離 (km) 的輔助函式 ---
get_pairwise_km <- function(sf_obj, label) {
  # 計算大圓距離矩陣 (公尺)
  dist_matrix <- st_distance(sf_obj)
  # 提取下三角矩陣 (排除對角線的 0 與對稱重複)，轉為公里 (km)
  dists_km <- as.numeric(dist_matrix[lower.tri(dist_matrix)]) / 1000
  data.frame(
    distance_km = dists_km,
    dataset = label
  )
}

# --- 2. 準備三組資料集的 Unique Points ---
# A. Moz 30 Clusters (30 個群聚中心代表點)
pts_moz30 <- moz_clusters_30_sf

# B. COMSA 116 Clusters (116 個獨立群聚聚合中心點)
pts_comsa116 <- comsa_sf

# C. Base AC (4,874 個獨立點位)
pts_baseac <- base_ac_sf

# --- 3. 計算成對距離並抽樣 Base AC (加速繪圖) ---
dist_moz30    <- get_pairwise_km(pts_moz30,    label = paste0("Moz Clusters (N = ", nrow(pts_moz30), ")"))
dist_comsa116 <- get_pairwise_km(pts_comsa116, label = paste0("COMSA Clusters (N = ", nrow(pts_comsa116), ")"))

# Base AC 有約 1,187 萬對，隨機抽樣 30 萬對保留完整分佈形態並優化渲染效能
dist_baseac_raw <- get_pairwise_km(pts_baseac, label = paste0("Base AC Points (N = ", nrow(pts_baseac), ")"))
set.seed(42)
dist_baseac <- if (nrow(dist_baseac_raw) > 300000) {
  dist_baseac_raw[sample(nrow(dist_baseac_raw), 300000), ]
} else {
  dist_baseac_raw
}

# --- 4. 合併資料並設定因子順序（由上到下：N=30 -> N=116 -> N=4874） ---
all_unique_dists <- bind_rows(dist_moz30, dist_comsa116, dist_baseac)
all_unique_dists$dataset <- factor(
  all_unique_dists$dataset,
  levels = c(
    paste0("Moz Clusters (N = ", nrow(pts_moz30), ")"),
    paste0("COMSA Clusters (N = ", nrow(pts_comsa116), ")"),
    paste0("Base AC Points (N = ", nrow(pts_baseac), ")")
  )
)

# --- 5. 繪製由上到下垂直堆疊的核密度分佈圖 ---
p_dist_vertical <- ggplot(all_unique_dists, aes(x = distance_km, fill = dataset, color = dataset)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  # ncol = 1 強制由上到下垂直排列；scales = "free_y" 保證不同樣本數下的密度峰值清晰可見
  facet_wrap(~ dataset, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("#2CA25F", "#E34A33", "#2B8CBE")) +
  scale_color_manual(values = c("#2CA25F", "#E34A33", "#2B8CBE")) +
  labs(
    title = "Pairwise Geodesic Distance Distribution Across Datasets",
    subtitle = "Calculated from unique spatial coordinates / centroids in Zambézia",
    x = "Pairwise Distance (km)",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray92", color = NA),
    panel.spacing = unit(1.2, "lines"),
    axis.title = element_text(face = "bold")
  )

print(p_dist_vertical)
