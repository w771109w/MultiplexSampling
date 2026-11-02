# ==============================================================================
# 專案：Zambézia 獨立病原體血清盛行率空間分佈圖 (紅藍雙變量 + 百分比 Cutoffs 標註)
# 輸入檔案：
#   - 01_data/processed/base_ac_processed.rds
# 輸出目錄：
#   - 03_results/maps_by_pathogen/
# 說明：
#   - 雙抗原病原 (7 隻)：自定義 9 色藍紅混色色盤 (X=藍, Y=紅, 雙高=深紫, 雙低=灰白)
#   - 單抗原病原 (3 隻)：Magma 連續漸層色階
#   - 地圖維持純淨，不常駐顯示縣名文字標籤
# ==============================================================================

library(here)
library(sf)
library(dplyr)
library(ggplot2)
library(cowplot)
library(scales)

# 確保已安裝並載入雙變量色階套件 biscale
if (!requireNamespace("biscale", quietly = TRUE)) {
  install.packages("biscale", repos = "https://cloud.r-project.org")
}
library(biscale)

# 建立獨立輸出資料夾
out_dir <- here("03_results", "maps_by_pathogen")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 關閉 s2 平面運算以防多邊形幾何自交報錯
sf_use_s2(FALSE)

# ------------------------------------------------------------------------------
# 步驟 1：載入 22 縣邊界多邊形與 30 個血清 Cluster
# ------------------------------------------------------------------------------
cat(">>> [1/3] 讀取 Zambézia 22 縣邊界多邊形與 30 個血清點位資料...\n")

# 1. 讀取主資料集並篩選 30 個 Sero 點位 (COMSA == 2)
base_ac_processed <- readRDS(here("01_data", "processed", "base_ac_processed.rds")) %>%
  st_transform(4326)

pts_sero30 <- base_ac_processed %>%
  filter(COMSA == 2)

cat("  - 30 個血清 Cluster 點位確認完畢，筆數:", nrow(pts_sero30), "\n")

# 2. 構建現行 22 縣邊界 (利用 GADM level 3 精確分離新設立縣)
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
  st_transform(4326)

cat("  - 22 縣多邊形幾何確認完畢，轄下 Distritos 總數:", nrow(zambezia_distritos), "\n")


# ------------------------------------------------------------------------------
# 步驟 2：核心通用繪圖函數 (自定義 9 色藍紅矩陣 + 百分比 Cutoffs 標註)
# ------------------------------------------------------------------------------
cat("\n>>> [2/3] 定義病原層級繪圖模組 (含自定義藍紅色盤與 Cutoffs 計算)...\n")

# 定義標準 3x3 藍紅自定義調色盤 (自下而上、由左至右對應 bi_class 1-1 到 3-3)
custom_blue_red <- c(
  "1-1" = "#E8E8E8", # 雙低: 淺灰白
  "2-1" = "#A4BED5", # X中: 淡灰藍
  "3-1" = "#1F77B4", # X高: 正深藍
  "1-2" = "#E8A3A8", # Y中: 淡粉紅
  "2-2" = "#A887A4", # 雙中: 灰紫
  "3-2" = "#62558E", # X高Y中: 藍紫
  "1-3" = "#D62728", # Y高: 正鮮紅
  "2-3" = "#9E3D6E", # Y高X中: 紅紫
  "3-3" = "#4A148C"  # 雙高: 深紫黑 (熱點)
)

# --- 2.1 雙抗原病原繪圖函數 (Bivariate Map with Cutoffs) ---
plot_pathogen_bivariate_cutoffs <- function(pts_sf, district_sf, 
                                            x_col, y_col, 
                                            x_label, y_label,
                                            pathogen_name, 
                                            file_prefix) {
  
  # 1. 預處理欄位名稱以規避 NSE 與向量檢查錯誤
  pts_prep <- pts_sf %>%
    mutate(
      x_val = .data[[x_col]],
      y_val = .data[[y_col]]
    )
  
  # 2. 精確計算三分位數 Cutoffs (33.3% 與 66.7%)
  x_cuts <- quantile(pts_prep$x_val, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  y_cuts <- quantile(pts_prep$y_val, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  
  fmt_pct <- function(val) paste0(round(val * 100, 1), "%")
  
  x_str <- paste0(
    x_label, " [X/藍]: ",
    "Low (<", fmt_pct(x_cuts[2]), ") | ",
    "Mid (", fmt_pct(x_cuts[2]), "–", fmt_pct(x_cuts[3]), ") | ",
    "High (≥", fmt_pct(x_cuts[3]), ")"
  )
  
  y_str <- paste0(
    y_label, " [Y/紅]: ",
    "Low (<", fmt_pct(y_cuts[2]), ") | ",
    "Mid (", fmt_pct(y_cuts[2]), "–", fmt_pct(y_cuts[3]), ") | ",
    "High (≥", fmt_pct(y_cuts[3]), ")"
  )
  
  # 3. 執行 3x3 雙變量分類
  pts_bi <- bi_class(
    pts_prep,
    x = x_val,
    y = y_val,
    style = "quantile",
    dim = 3
  )
  
  # 4. 繪製主地圖 (不顯示縣名文字標籤，傳入自定義 custom_blue_red 色盤)
  p_map <- ggplot() +
    geom_sf(
      data = district_sf, 
      fill = "#F8F9FA", 
      color = "#90A4AE", 
      linewidth = 0.4
    ) +
    geom_sf(
      data = pts_bi,
      aes(fill = bi_class),
      shape = 21,
      size = 4.8,
      stroke = 0.75,
      color = "#263238"
    ) +
    bi_scale_fill(pal = custom_blue_red, dim = 3) +
    labs(
      title = paste0(pathogen_name, " Dual-Antigen Bivariate Map"),
      subtitle = paste0(x_str, "\n", y_str),
      caption = "Source: Mozambique Age Serosurvey & COMSA Frame (Zambézia 22 Distritos, N=30 Clusters)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold.italic", size = 13.5, hjust = 0),
      plot.subtitle = element_text(color = "#37474F", size = 8.8, lineheight = 1.25, margin = margin(b = 8)),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.position = "none" # 主地圖不放預設圖例，由 cowplot 疊加 2D 圖例
    )
  
  # 5. 建立 3x3 二維圖例 (自訂標籤為 X=藍, Y=紅)
  p_legend <- bi_legend(
    pal = custom_blue_red,
    dim = 3,
    xlab = paste0("Higher ", x_label, " (Blue)"),
    ylab = paste0("Higher ", y_label, " (Red)"),
    size = 7
  )
  
  # 6. 拼裝地圖與二維圖例 (將圖例置於右下角外海空白處)
  p_combined <- ggdraw() +
    draw_plot(p_map, 0, 0, 1, 1) +
    draw_plot(p_legend, 0.69, 0.11, 0.27, 0.27)
  
  # 輸出高解析檔案
  out_png <- file.path(out_dir, paste0("map_bivariate_", file_prefix, ".png"))
  out_pdf <- file.path(out_dir, paste0("map_bivariate_", file_prefix, ".pdf"))
  
  ggsave(out_png, p_combined, width = 9.8, height = 7.5, dpi = 300)
  ggsave(out_pdf, p_combined, width = 9.8, height = 7.5)
  
  cat("  ✓ [雙抗原病原] ", pathogen_name, "\n")
  cat("     -", x_str, "\n")
  cat("     -", y_str, "\n")
}

# --- 2.2 單抗原病原繪圖函數 (Continuous Gradient Map) ---
plot_pathogen_single <- function(pts_sf, district_sf,
                                 antigen_col, antigen_label,
                                 pathogen_name,
                                 file_prefix) {
  
  # 繪製主地圖 (不顯示縣名文字標籤)
  p_single <- ggplot() +
    geom_sf(
      data = district_sf, 
      fill = "#F8F9FA", 
      color = "#90A4AE", 
      linewidth = 0.4
    ) +
    geom_sf(
      data = pts_sf,
      aes(color = .data[[antigen_col]], fill = .data[[antigen_col]]),
      shape = 21,
      size = 4.8,
      stroke = 0.75
    ) +
    scale_color_viridis_c(
      option = "magma",
      direction = -1,
      name = "Seroprevalence",
      labels = percent_format(accuracy = 1)
    ) +
    scale_fill_viridis_c(
      option = "magma",
      direction = -1,
      name = "Seroprevalence",
      labels = percent_format(accuracy = 1)
    ) +
    labs(
      title = paste0(pathogen_name, " (", antigen_label, ")"),
      subtitle = "Single Antigen Seroprevalence (Zambézia Province, 22 Distritos, N = 30 Sero Clusters)",
      caption = "Source: Mozambique Age Serosurvey & COMSA Sampling Frame"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold.italic", size = 13.5, hjust = 0),
      plot.subtitle = element_text(color = "#455A64", size = 9.5, margin = margin(b = 8)),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.position = "right",
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.4, "cm")
    )
  
  out_png <- file.path(out_dir, paste0("map_single_", file_prefix, ".png"))
  out_pdf <- file.path(out_dir, paste0("map_single_", file_prefix, ".pdf"))
  
  ggsave(out_png, p_single, width = 9, height = 7.5, dpi = 300)
  ggsave(out_pdf, p_single, width = 9, height = 7.5)
  cat("  ✓ [單抗原病原] ", pathogen_name, "已完成匯出。\n")
}


# ------------------------------------------------------------------------------
# 步驟 3：依生物學病原體 (Pathogen Species) 逐一執行並匯出
# ------------------------------------------------------------------------------
cat("\n>>> [3/3] 開始繪製 7 隻雙抗原病原體 + 3 隻單抗原病原地圖...\n")

# --- A. 雙抗原病原體 (共 7 隻，紅藍色階 + 標註低中高 Cutoffs) ---

# 1. Brugia malayi (馬來絲蟲)
plot_pathogen_bivariate_cutoffs(
  pts_sf = pts_sero30, district_sf = zambezia_distritos,
  x_col = "seroprev_bm14", y_col = "seroprev_bm33",
  x_label = "Bm14", y_label = "Bm33",
  pathogen_name = "Brugia malayi",
  file_prefix = "01_brugia_malayi"
)

# 2. Chlamydia trachomatis (砂眼披衣菌)
plot_pathogen_bivariate_cutoffs(
  pts_sf = pts_sero30, district_sf = zambezia_distritos,
  x_col = "seroprev_pgp3", y_col = "seroprev_ct694",
  x_label = "Pgp3", y_label = "Ct694",
  pathogen_name = "Chlamydia trachomatis",
  file_prefix = "02_chlamydia_trachomatis"
)

# 3. Schistosoma mansoni (曼森血吸蟲)
plot_pathogen_bivariate_cutoffs(
  pts_sf = pts_sero30, district_sf = zambezia_distritos,
  x_col = "seroprev_sm25", y_col = "seroprev_sea",
  x_label = "Sm25", y_label = "SEA",
  pathogen_name = "Schistosoma mansoni",
  file_prefix = "03_schistosoma_mansoni"
)

# 4. Cryptosporidium spp. (隱孢子蟲)
plot_pathogen_bivariate_cutoffs(
  pts_sf = pts_sero30, district_sf = zambezia_distritos,
  x_col = "seroprev_cp17", y_col = "seroprev_cp23",
  x_label = "Cp17", y_label = "Cp23",
  pathogen_name = "Cryptosporidium spp.",
  file_prefix = "04_cryptosporidium_spp"
)

# 5. Giardia lamblia (梨形鞭毛蟲)
plot_pathogen_bivariate_cutoffs(
  pts_sf = pts_sero30, district_sf = zambezia_distritos,
  x_col = "seroprev_vsp3", y_col = "seroprev_vsp5",
  x_label = "Vsp3", y_label = "Vsp5",
  pathogen_name = "Giardia lamblia",
  file_prefix = "05_giardia_lamblia"
)

# 6. Taenia solium (豬肉絛蟲 / 囊蟲)
plot_pathogen_bivariate_cutoffs(
  pts_sf = pts_sero30, district_sf = zambezia_distritos,
  x_col = "seroprev_t24h", y_col = "seroprev_es33",
  x_label = "T24h", y_label = "ES33",
  pathogen_name = "Taenia solium",
  file_prefix = "06_taenia_solium"
)

# 7. Treponema pallidum (梅毒 / 雅司密螺旋體)
plot_pathogen_bivariate_cutoffs(
  pts_sf = pts_sero30, district_sf = zambezia_distritos,
  x_col = "seroprev_rp17", y_col = "seroprev_tmpa",
  x_label = "Rp17", y_label = "TmpA",
  pathogen_name = "Treponema pallidum",
  file_prefix = "07_treponema_pallidum"
)


# --- B. 單一抗原病原體 (共 3 隻，連續漸層色階) ---

# 8. Wuchereria bancrofti (班氏絲蟲)
plot_pathogen_single(
  pts_sf = pts_sero30, district_sf = zambezia_distritos,
  antigen_col = "seroprev_wb123", antigen_label = "Wb123",
  pathogen_name = "Wuchereria bancrofti",
  file_prefix = "08_wuchereria_bancrofti"
)

# 9. Strongyloides stercoralis (糞線蟲)
plot_pathogen_single(
  pts_sf = pts_sero30, district_sf = zambezia_distritos,
  antigen_col = "seroprev_nie", antigen_label = "NIE",
  pathogen_name = "Strongyloides stercoralis",
  file_prefix = "09_strongyloides_stercoralis"
)

# 10. Onchocerca volvulus (旋盤尾絲蟲 / 河盲)
plot_pathogen_single(
  pts_sf = pts_sero30, district_sf = zambezia_distritos,
  antigen_col = "seroprev_ov16", antigen_label = "Ov16",
  pathogen_name = "Onchocerca volvulus",
  file_prefix = "10_onchocerca_volvulus"
)

# 還原 s2 球面幾何引擎設定
sf_use_s2(TRUE)

cat("\n==============================================================================\n")
cat("全數完成！10 隻病原地圖 (PNG 與 PDF) 已成功匯出至：\n")
cat(out_dir, "\n")
cat("==============================================================================\n")