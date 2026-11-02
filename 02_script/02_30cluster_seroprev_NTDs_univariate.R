# ==============================================================================
# 專案：Zambézia 15 種抗原獨立連續色階空間地圖 (單抗原獨立出圖)
# 輸入檔案：
#   - 01_data/processed/base_ac_processed.rds
# 輸出目錄：
#   - 03_results/maps_single_antigens/ (每種抗原一張 PNG 與 PDF)
# ==============================================================================

library(here)
library(sf)
library(dplyr)
library(ggplot2)
library(scales)

# 建立獨立專屬輸出資料夾
out_single_dir <- here("03_results", "maps_single_antigens")
dir.create(out_single_dir, recursive = TRUE, showWarnings = FALSE)

# 關閉 s2 平面運算以防多邊形自交報錯
sf_use_s2(FALSE)

# ------------------------------------------------------------------------------
# 步驟 1：載入 22 縣邊界多邊形與 30 個血清 Cluster
# ------------------------------------------------------------------------------
cat(">>> [1/3] 讀取 Zambézia 22 縣行政區邊界與 30 個血清點位資料...\n")

# 1. 讀取主資料集並篩選 30 個 Sero 代表點 (COMSA == 2)
base_ac_processed <- readRDS(here("01_data", "processed", "base_ac_processed.rds")) %>%
  st_transform(4326)

pts_sero30 <- base_ac_processed %>%
  filter(COMSA == 2)

cat("  - 30 個血清 Cluster 點位確認完畢，筆數:", nrow(pts_sero30), "\n")

# 2. 構建現行 22 縣行政區多邊形 (利用 GADM level 3 精確分離新設立縣)
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
# 步驟 2：定義 15 種抗原的病原體分類與詳細標籤字典
# ------------------------------------------------------------------------------
cat("\n>>> [2/3] 配置 15 種抗原的中英文元數據清單...\n")

antigen_meta <- list(
  # 淋巴絲蟲病 (Lymphatic Filariasis)
  list(col = "seroprev_wb123", name = "Wb123", pathogen = "Wuchereria bancrofti (Lymphatic Filariasis)", note = "Recent / Active Transmission Marker"),
  list(col = "seroprev_bm14",  name = "Bm14",  pathogen = "Brugia malayi / W. bancrofti",             note = "Historical Cumulative Exposure Marker"),
  list(col = "seroprev_bm33",  name = "Bm33",  pathogen = "Brugia malayi / W. bancrofti",             note = "High-Sensitivity Filariasis Exposure"),
  
  # 砂眼 (Trachoma)
  list(col = "seroprev_pgp3",  name = "Pgp3",  pathogen = "Chlamydia trachomatis (Trachoma)",         note = "Major Outer Membrane / Active Exposure"),
  list(col = "seroprev_ct694", name = "Ct694", pathogen = "Chlamydia trachomatis (Trachoma)",         note = "Secreted Early Pathogenesis Antigen"),
  
  # 曼森血吸蟲病 (Schistosomiasis)
  list(col = "seroprev_sm25",  name = "Sm25",  pathogen = "Schistosoma mansoni (Schistosomiasis)",    note = "Adult Worm Surface Glycoprotein"),
  list(col = "seroprev_sea",   name = "SEA",   pathogen = "Schistosoma mansoni (Schistosomiasis)",    note = "Soluble Egg Antigen (Intensity/Fecundity)"),
  
  # 水傳播原蟲 (Waterborne Protozoa / WASH)
  list(col = "seroprev_cp17",  name = "Cp17",  pathogen = "Cryptosporidium spp.",                     note = "Sporozoite Surface Antigen (WASH)"),
  list(col = "seroprev_cp23",  name = "Cp23",  pathogen = "Cryptosporidium spp.",                     note = "Oocyst Wall & Sporozoite Antigen (WASH)"),
  list(col = "seroprev_vsp3",  name = "Vsp3",  pathogen = "Giardia lamblia",                          note = "Variant-Specific Surface Protein 3 (WASH)"),
  list(col = "seroprev_vsp5",  name = "Vsp5",  pathogen = "Giardia lamblia",                          note = "Variant-Specific Surface Protein 5 (WASH)"),
  
  # 絛蟲 / 囊蟲症 (Taeniasis / Cysticercosis)
  list(col = "seroprev_t24h",  name = "T24h",  pathogen = "Taenia solium (Cysticercosis)",            note = "Diagnostic Tetraspanin for Neurocysticercosis"),
  list(col = "seroprev_es33",  name = "ES33",  pathogen = "Taenia solium (Cysticercosis)",            note = "Excretory-Secretory Antigen"),
  
  # 梅毒 / 雅司密螺旋體 (Treponema)
  list(col = "seroprev_rp17",  name = "Rp17",  pathogen = "Treponema pallidum (Yaws / Syphilis)",     note = "Membrane Lipoprotein Antigen"),
  list(col = "seroprev_tmpa",  name = "TmpA",  pathogen = "Treponema pallidum (Yaws / Syphilis)",     note = "Treponemal Membrane Protein A"),
  
  # 糞線蟲症 (Strongyloidiasis)
  list(col = "seroprev_nie",   name = "NIE",   pathogen = "Strongyloides stercoralis",                note = "Larval Recombinant Antigen"),
  
  # 盤尾絲蟲病 / 河盲症 (Onchocerciasis)
  list(col = "seroprev_ov16",  name = "Ov16",  pathogen = "Onchocerca volvulus (River Blindness)",    note = "IgG4 Recombinant Antigen (Elimination Benchmark)")
)


# ------------------------------------------------------------------------------
# 步驟 3：批次產出 15+ 張單抗原獨立地圖至新資料夾 (無文字純淨圖面)
# ------------------------------------------------------------------------------
cat("\n>>> [3/3] 開始批次繪製 15 種抗原獨立空間地圖並儲存...\n")

for (i in seq_along(antigen_meta)) {
  item <- antigen_meta[[i]]
  col_name <- item$col
  
  # 檢查該欄位是否存在於資料集中
  if (!col_name %in% names(pts_sero30)) {
    warning("找不到欄位: ", col_name, "，略過本圖。")
    next
  }
  
  # 計算該抗原在 30 點的統計分佈
  vals <- pts_sero30[[col_name]]
  val_min <- min(vals, na.rm = TRUE)
  val_med <- median(vals, na.rm = TRUE)
  val_max <- max(vals, na.rm = TRUE)
  
  subtitle_text <- paste0(
    "Target: ", item$name, " (", item$note, ")\n",
    "Across 30 Clusters: Median = ", round(val_med * 100, 1), "% | Range = [",
    round(val_min * 100, 1), "% – ", round(val_max * 100, 1), "%]"
  )
  
  # 繪製主地圖 (不顯示縣名文字)
  p <- ggplot() +
    # 22 縣行政區多邊形底圖
    geom_sf(
      data = zambezia_distritos,
      fill = "#F8F9FA",
      color = "#90A4AE",
      linewidth = 0.4
    ) +
    # 30 個 Cluster 點位 (連續色階漸層)
    geom_sf(
      data = pts_sero30,
      aes(color = .data[[col_name]], fill = .data[[col_name]]),
      shape = 21,
      size = 4.8,
      stroke = 0.75
    ) +
    scale_color_viridis_c(
      option = "magma",
      direction = -1,
      name = "Seroprevalence",
      labels = percent_format(accuracy = 1),
      limits = c(0, max(0.01, val_max))
    ) +
    scale_fill_viridis_c(
      option = "magma",
      direction = -1,
      name = "Seroprevalence",
      labels = percent_format(accuracy = 1),
      limits = c(0, max(0.01, val_max))
    ) +
    labs(
      title = paste0(item$pathogen, " — ", item$name),
      subtitle = subtitle_text,
      caption = "Source: Mozambique Age Serosurvey & COMSA Frame (Zambézia 22 Distritos, N=30 Clusters)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = element_text(color = "#37474F", size = 9, lineheight = 1.25, margin = margin(b = 10)),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.position = "right",
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.4, "cm")
    )
  
  # 輸出高解析 PNG 與 PDF
  file_slug <- sprintf("%02d_%s", i, sub("^seroprev_", "", col_name))
  out_png <- file.path(out_single_dir, paste0("map_single_", file_slug, ".png"))
  out_pdf <- file.path(out_single_dir, paste0("map_single_", file_slug, ".pdf"))
  
  ggsave(out_png, p, width = 9.2, height = 7.5, dpi = 300)
  ggsave(out_pdf, p, width = 9.2, height = 7.5)
  
  cat(sprintf("  [%02d/17] ✓ 已儲存: %s (%s)\n", i, item$name, item$pathogen))
}

# 還原 s2 球面幾何引擎設定
sf_use_s2(TRUE)

cat("\n==============================================================================\n")
cat("全數完成！單抗原獨立地圖已全數匯出至專屬目錄：\n")
cat(out_single_dir, "\n")
cat("==============================================================================\n")