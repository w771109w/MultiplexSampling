# ==============================================================================
# 專案：Zambézia 省 22 縣 17 種抗原血清盛行率基準估計 (Wilson 95% CI)
#       與空間地圖 + 森林圖 (Map + Forest Plot) 複合面板多頁 Atlas
# 輸入：
#   - 01_data/processed/base_ac_processed.rds
#   - 01_data/raw/moz_age.rds
# 輸出：
#   - 03_results/table_district_all_antigens_wilson_ci.csv
#   - 03_results/antigen_panels/ (個別 PNG 複合面板)
#   - 03_results/all_antigens_seroprevalence_atlas.pdf (整合 17 頁多頁 PDF)
# ==============================================================================

library(here)
library(sf)
library(dplyr)
library(tidyr)
library(scales)
library(ggplot2)
library(patchwork)
library(ggspatial)

# 關閉 s2 平面運算以確保多邊形幾何拓撲穩定
sf_use_s2(FALSE)

dir.create(here("03_results"), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 步驟 1：建立現行 22 縣邊界多邊形 (GADM Level 3 解析重整)
# ==============================================================================
cat(">>> [1/5] 下載並重構 Zambézia 現行 22 縣行政邊界...\n")

if (!requireNamespace("geodata", quietly = TRUE)) install.packages("geodata")
moz_adm3_raw <- geodata::gadm(country = "MOZ", level = 3, path = tempdir()) %>%
  st_as_sf()

zambezia_distritos_22 <- moz_adm3_raw %>%
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
  rename(district_gadm = district_22) %>%
  mutate(district_gadm = trimws(district_gadm))

# ==============================================================================
# 步驟 2：直接透過「空間交集」將 30 個 Moz Clusters 綁定到 22 縣 (核心修復點)
# ==============================================================================
cat("\n>>> [2/5] 讀取 moz_age.rds 並以空間經緯度精確對齊 22 縣...\n")

moz_age <- readRDS(here("01_data", "raw", "moz_age.rds"))

# 2.1 以每個 csa_code 的家戶經緯度平均，建立 30 個群聚的空間點位
moz_clusters_sf <- moz_age %>%
  group_by(csa_code) %>%
  summarize(
    mean_lon = mean(lon, na.rm = TRUE),
    mean_lat = mean(lat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  st_as_sf(coords = c("mean_lon", "mean_lat"), crs = 4326)

# 2.2 與 22 縣多邊形進行空間交集 (Spatial Join)
moz_clusters_aligned <- st_join(moz_clusters_sf, zambezia_distritos_22, join = st_intersects)

# 邊界外點位微調
na_pts <- which(is.na(moz_clusters_aligned$district_gadm))
if (length(na_pts) > 0) {
  nearest_idx <- st_nearest_feature(moz_clusters_aligned[na_pts, ], zambezia_distritos_22)
  moz_clusters_aligned$district_gadm[na_pts] <- zambezia_distritos_22$district_gadm[nearest_idx]
}

# 建立 100% 絕對正確的對照表：csa_code -> district_gadm
csa_to_district <- moz_clusters_aligned %>%
  st_drop_geometry() %>%
  select(csa_code, district_gadm) %>%
  mutate(csa_code = as.character(csa_code))

cat("  - 30 個 Sero 群聚分佈縣份檢驗:\n")
print(table(csa_to_district$district_gadm))

# 2.3 將 1,292 個體資料對齊行政區
moz_age_aligned <- moz_age %>%
  mutate(csa_code = as.character(csa_code)) %>%
  left_join(csa_to_district, by = "csa_code")

cat("  ✓ 1,292 個體成功對齊行政區數筆:", sum(!is.na(moz_age_aligned$district_gadm)), "/ 1292\n")

# ==============================================================================
# 步驟 3：讀取 1,292 個人資料，批次計算全部抗原之 Wilson Score 95% CI
# ==============================================================================
cat("\n>>> [3/5] 批次計算全部抗原之 Wilson Score Interval...\n")

pos_cols <- grep("^pos_", names(moz_age), value = TRUE)
z_val <- qnorm(0.975)

# 行政區層級計算 (12 縣)
district_counts <- moz_age_aligned %>%
  filter(!is.na(district_gadm)) %>%
  pivot_longer(
    cols = all_of(pos_cols),
    names_to = "antigen",
    names_prefix = "pos_",
    values_to = "pos_status"
  ) %>%
  group_by(district_gadm, antigen) %>%
  summarize(
    sample_n = sum(!is.na(pos_status)),
    pos_n    = sum(pos_status == 1, na.rm = TRUE),
    .groups  = "drop"
  )

district_all_antigens_long <- district_counts %>%
  mutate(
    p = if_else(sample_n > 0, pos_n / sample_n, NA_real_),
    denom = 1 + (z_val^2) / sample_n,
    center = (p + (z_val^2) / (2 * sample_n)) / denom,
    margin = (z_val * sqrt((p * (1 - p) / sample_n) + (z_val^2) / (4 * sample_n^2))) / denom,
    mean  = p,
    lower = pmax(0, center - margin),
    upper = pmin(1, center + margin),
    seroprev_formatted = sprintf("%.1f%% (%.1f%%–%.1f%%)", mean * 100, lower * 100, upper * 100)
  )

# 全省總計列
provincial_counts <- moz_age_aligned %>%
  pivot_longer(
    cols = all_of(pos_cols),
    names_to = "antigen",
    names_prefix = "pos_",
    values_to = "pos_status"
  ) %>%
  group_by(antigen) %>%
  summarize(
    district_gadm = "Total (Zambézia)",
    sample_n      = sum(!is.na(pos_status)),
    pos_n         = sum(pos_status == 1, na.rm = TRUE),
    .groups       = "drop"
  )

provincial_all_antigens <- provincial_counts %>%
  mutate(
    p = if_else(sample_n > 0, pos_n / sample_n, NA_real_),
    denom = 1 + (z_val^2) / sample_n,
    center = (p + (z_val^2) / (2 * sample_n)) / denom,
    margin = (z_val * sqrt((p * (1 - p) / sample_n) + (z_val^2) / (4 * sample_n^2))) / denom,
    mean  = p,
    lower = pmax(0, center - margin),
    upper = pmin(1, center + margin),
    seroprev_formatted = sprintf("%.1f%% (%.1f%%–%.1f%%)", mean * 100, lower * 100, upper * 100)
  )

all_antigens_long_final <- bind_rows(district_all_antigens_long, provincial_all_antigens)

# 匯出寬表格
district_all_antigens_wide <- all_antigens_long_final %>%
  select(district_gadm, antigen, seroprev_formatted) %>%
  pivot_wider(names_from = antigen, values_from = seroprev_formatted)

write.csv(
  district_all_antigens_wide,
  here("03_results", "table_district_all_antigens_wilson_ci.csv"),
  row.names = FALSE
)

# ==============================================================================
# 步驟 4：批次產出方案 3 複合面板 (Panel A: 地圖 + Panel B: 森林圖)
# ==============================================================================
cat("\n>>> [4/5] 批次繪製地圖 + 森林圖複合面板並輸出至 PDF Atlas...\n")

# 17 種抗原對應之病原體與疾病字典表
antigen_metadata <- tibble::tribble(
  ~antigen, ~pathogen,                    ~disease,
  "bm14",   "Brugia malayi",       "Lymphatic Filariasis",
  "bm33",   "Brugia malayi",       "Lymphatic Filariasis",
  "wb123",  "Wuchereria bancrofti",       "Lymphatic Filariasis",
  "sm25",   "Schistosoma mansoni",        "Schistosomiasis",
  "sea",    "Schistosoma mansoni",        "Schistosomiasis",
  "pgp3",   "Chlamydia trachomatis",      "Trachoma",
  "ct694",  "Chlamydia trachomatis",      "Trachoma",
  "t24h",   "Taenia solium",              "Cysticercosis / Taeniasis",
  "es33",   "Taenia solium",              "Cysticercosis / Taeniasis",
  "ov16",   "Onchocerca volvulus",        "Onchocerciasis (River Blindness)",
  "tmpa",   "Treponema pallidum",         "Syphilis / Yaws",
  "rp17",   "Treponema pallidum",         "Syphilis / Yaws",
  "nie",    "Strongyloides stercoralis",  "Strongyloidiasis",
  "vsp5",   "Giardia duodenalis",         "Giardiasis",
  "vsp3",   "Giardia duodenalis",         "Giardiasis",
  "cp23",   "Cryptosporidium parvum",     "Cryptosporidiosis",
  "cp17",   "Cryptosporidium parvum",     "Cryptosporidiosis"
)

maps_png_dir <- here("03_results", "antigen_panels")
dir.create(maps_png_dir, recursive = TRUE, showWarnings = FALSE)

antigen_list <- sub("^pos_", "", pos_cols)

pdf_path <- here("03_results", "all_antigens_seroprevalence_atlas.pdf")
pdf(pdf_path, width = 14, height = 7.5, onefile = TRUE)

for (ag in antigen_list) {
  cat("  - 正在產製複合面板:", ag, "...\n")
  
  meta_info <- antigen_metadata %>% filter(antigen == ag)
  if (nrow(meta_info) > 0) {
    pathogen_name <- meta_info$pathogen[1]
    disease_name  <- meta_info$disease[1]
  } else {
    pathogen_name <- "Multiplex Antigen"
    disease_name  <- "Serological Marker"
  }
  
  # 1. 取得該抗原各縣統計值
  ag_district_stats <- district_all_antigens_long %>%
    filter(antigen == ag) %>%
    select(district_gadm, mean, lower, upper, seroprev_formatted, sample_n)
  
  # 2. 取得全省統計值
  ag_prov_stats <- provincial_all_antigens %>% filter(antigen == ag)
  prov_mean  <- ag_prov_stats$mean
  prov_lower <- ag_prov_stats$lower
  prov_upper <- ag_prov_stats$upper
  
  prov_text <- sprintf("Provincial Mean: %.1f%% (95%% CI: %.1f%%–%.1f%%)",
                       prov_mean * 100, prov_lower * 100, prov_upper * 100)
  
  # 3. 地圖資料連結
  ag_map_sf <- zambezia_distritos_22 %>%
    left_join(ag_district_stats, by = "district_gadm") %>%
    mutate(mean_pct = mean * 100)
  
  # 色階上限自適應 (動態調整，保證各縣能充分展現顏色漸層)
  max_val <- max(ag_map_sf$mean_pct, na.rm = TRUE)
  upper_limit <- if (is.finite(max_val) && max_val > 0) ceiling(max_val * 1.1) else 20
  upper_limit <- max(upper_limit, ceiling(prov_mean * 100 * 1.15))
  
  # 4. 繪製 Panel A：空間地圖 (Choropleth)
  p_map <- ggplot() +
    geom_sf(
      data = ag_map_sf,
      aes(fill = mean_pct),
      color = "#424242",
      linewidth = 0.45
    ) +
    geom_sf(
      data = moz_clusters_aligned,
      aes(shape = "Serosurvey Site (N=30)"),
      fill = "#E53935",
      color = "black",
      size = 2.4,
      stroke = 0.45
    ) +
    scale_fill_viridis_c(
      name = "Seroprevalence (%)",
      option = "plasma",
      na.value = "#E0E0E0", # 無抽樣的 10 個縣填灰色
      limits = c(0, upper_limit),
      oob = scales::squish,
      labels = function(x) paste0(x, "%"),
      guide = guide_colorbar(
        barwidth = 0.8,
        barheight = 8.5,
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    scale_shape_manual(name = "", values = c("Serosurvey Site (N=30)" = 21)) +
    annotation_scale(location = "bl", width_hint = 0.22, text_size = 8) +
    annotation_north_arrow(
      location = "tr",
      which_north = "true",
      height = unit(1, "cm"),
      width = unit(1, "cm"),
      style = north_arrow_fancy_orienteering
    ) +
    labs(
      title = "A  Spatial Seroprevalence Distribution",
      subtitle = "Grey regions denote 10 districts without serosurvey sampling",
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(color = "grey35", size = 9),
      panel.background = element_rect(fill = "#F4F7F9", color = NA),
      panel.grid.major = element_line(color = "white", linewidth = 0.35),
      legend.position = "right"
    )
  
  # 5. 繪製 Panel B：森林圖 (依 Mean 由低到高排列 factor level，使畫布呈現由上到下「高到低」)
  forest_df <- ag_district_stats %>%
    filter(!is.na(mean) & sample_n > 0) %>%
    # arrange(mean) 配合 factor levels，會在 ggplot Y 軸上呈現「由上至下：高 -> 低」
    arrange(mean) %>%
    mutate(
      district_ordered = factor(district_gadm, levels = district_gadm),
      pct_label = sprintf("%.1f%% (N=%d)", mean * 100, sample_n)
    )
  
  p_forest <- ggplot(forest_df, aes(y = district_ordered, x = mean)) +
    # 全省平均虛線
    geom_vline(xintercept = prov_mean, linetype = "dashed", color = "#D50000", linewidth = 0.75) +
    # 95% CI 橫向誤差線
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25, color = "#283593", linewidth = 0.85) +
    # 點估計值 (依樣本數設大小)
    geom_point(aes(size = sample_n), color = "#0D47A1", fill = "#42A5F5", shape = 21) +
    # 點右側標示百分比與樣本數
    geom_text(aes(label = pct_label), vjust = -0.7, size = 3.1, color = "#1A237E", fontface = "bold") +
    scale_x_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, max(c(forest_df$upper, prov_upper, 0.15), na.rm = TRUE) * 1.15),
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_size_continuous(name = "Sample (N)", range = c(2.5, 5)) +
    labs(
      title = "B  District-Level Estimates & Wilson 95% CI",
      subtitle = paste0("Ranked from highest to lowest | Dashed red line = Provincial Mean (", round(prov_mean * 100, 1), "%)"),
      x = "Seroprevalence",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(color = "grey35", size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linetype = "dotted"),
      axis.text.y = element_text(face = "bold", color = "#212121", size = 10),
      legend.position = "right"
    )
  
  # 6. 整合複合面板 (Patchwork)
  p_combined <- (p_map | p_forest) +
    plot_layout(widths = c(1.15, 1)) +
    plot_annotation(
      title = paste0("Antigen: ", toupper(ag), " (Target: ", pathogen_name, " | ", disease_name, ")"),
      subtitle = paste0(prov_text, "\nCalculated via Wilson Score 95% CI (Independent Sample Assumption, N = 1,292)"),
      caption = "Notes: Grey regions indicate 10 unsampled districts. Error bars denote 95% Wilson score confidence intervals.",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, color = "#0D47A1"),
        plot.subtitle = element_text(size = 10, color = "#37474F", lineheight = 1.2),
        plot.caption = element_text(size = 8.5, color = "grey40")
      )
    )
  
  # 7. 輸出 PNG
  png_path <- file.path(maps_png_dir, paste0("panel_seroprev_", ag, ".png"))
  ggsave(
    filename = png_path,
    plot = p_combined,
    width = 14,
    height = 7.5,
    dpi = 300
  )
  
  # 8. 印入多頁 PDF
  print(p_combined)
}

dev.off()
sf_use_s2(TRUE)

cat("\n==============================================================================\n")
cat("全數修正完成！\n")
cat("1. 17 種抗原複合面板 PNG 已輸出至:\n   03_results/antigen_panels/\n")
cat("2. 整合之 17 頁多頁 PDF Atlas 已儲存至:\n  ", pdf_path, "\n")
cat("==============================================================================\n")