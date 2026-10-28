# ==============================================================================
# 專案：莫三比克 Zambézia 省 4,874 取樣點抽樣分層互動式地圖 (Street & Satellite)
# 資料來源：01_data/processed/base_ac_processed.rds
# 輸出：RStudio Viewer 互動式地圖 (包含向量底圖與高解析衛星空照圖切換)
# ==============================================================================

library(here)
library(sf)
library(dplyr)
library(mapview)
library(leafpop)

# ------------------------------------------------------------------------------
# 1. 讀取處理後的 RDS 資料並建立英文化抽樣階層欄位
# ------------------------------------------------------------------------------
cat(">>> 讀取 base_ac_processed.rds 並建立圖層屬性標籤...\n")

base_ac_processed <- readRDS(here("01_data", "processed", "base_ac_processed.rds")) %>%
  st_transform(4326)

# 依 COMSA 欄位數值建立嚴謹的英文學術標籤與因子順序
base_ac_map <- base_ac_processed %>%
  mutate(
    # 建立英文分類標籤
    comsa_tier = case_when(
      COMSA == 2 ~ "COMSA + Serosurvey (N = 30)",
      COMSA == 1 ~ "COMSA Only (N = 85)",
      COMSA == 0 ~ "Baseline Sampling Frame (N = 4,759)",
      TRUE       ~ "Unclassified"
    ),
    # 固定圖例排列順序 (由最核心採樣至母體底層)
    comsa_tier = factor(
      comsa_tier,
      levels = c(
        "COMSA + Serosurvey (N = 30)",
        "COMSA Only (N = 85)",
        "Baseline Sampling Frame (N = 4,759)"
      )
    ),
    # 滑鼠懸停即時提示文字 (Tooltip)
    tooltip_text = paste0(
      "Tier: ", comsa_tier,
      if_else(!is.na(seroprev), paste0(" | Seroprev: ", round(seroprev, 3)), ""),
      if_else(!is.na(comsa_cluster_id), paste0(" | Cluster ID: ", comsa_cluster_id), "")
    )
  )

# 印出檢查數量
cat("4,874 點位階層分佈確認：\n")
print(table(base_ac_map$comsa_tier))

# ------------------------------------------------------------------------------
# 2. 設定圖磚來源：向量地圖 (Street) 與真實衛星空照圖 (Satellite)
# ------------------------------------------------------------------------------
mapviewOptions(
  basemaps = c(
    "OpenStreetMap",         # 標準向量街道與聚落底圖 (General Street Map)
    "Esri.WorldImagery",     # 高解析度真實衛星空照圖 (Satellite Imagery)
    "CartoDB.Positron"       # 極簡灰白向量底圖 (突顯點位顏色分佈)
  ),
  fgb = TRUE                 # FlatGeobuf 啟用，加速 4,800+ 點位在瀏覽器中的縮放效能
)

# 定義各階層在深淺底圖上皆清晰可辨的專屬色盤
tier_colors <- c(
  "COMSA + Serosurvey (N = 30)"         = "#00E676", # 亮翠綠 (最核心血清調查點)
  "COMSA Only (N = 85)"                = "#FF9100", # 明亮暖橘 (常規 COMSA 抽樣點)
  "Baseline Sampling Frame (N = 4,759)"= "#2979FF"  # 電光藍 (背景未選入母體點)
)

# ------------------------------------------------------------------------------
# 3. 繪製互動地圖 (mapview)
# ------------------------------------------------------------------------------
# 為了避免 4,759 個背景點遮擋住核心的 30 個血清點與 85 個 COMSA 點，
# 將資料依 COMSA 數值升冪排序 (先畫 0，再畫 1，最後畫 2 疊在最頂層)
base_ac_map_ordered <- base_ac_map %>%
  arrange(COMSA)

# 找出要呈現在彈出視窗 (Popup) 的關鍵欄位
popup_cols <- intersect(
  c("COMSA", "comsa_tier", "comsa_cluster_id", "seroprev", "EA_CODE", "DIST_CODE", "LATITUDE", "LONGITUDE"),
  names(base_ac_map_ordered)
)

map_interactive <- mapview(
  base_ac_map_ordered,
  zcol = "comsa_tier",
  col.regions = tier_colors,
  color = "black",              # 黑色外圈輪廓，無論在灰白街景或深綠樹林衛星圖中皆能突顯邊界
  lwd = 0.5,
  cex = 3.5,                    # 點位大小適中
  alpha.regions = 0.85,         # 半透明避免完全遮蔽衛星地貌
  legend = TRUE,
  label = base_ac_map_ordered$tooltip_text,
  layer.name = "Sampling Tiers (N = 4,874)",
  popup = popupTable(base_ac_map_ordered, zcol = popup_cols)
)

# 在 RStudio Viewer 面板中渲染呈現
print(map_interactive)
cat("\n>>> 互動式地圖已成功渲染至 RStudio Viewer 面板！\n")