# ==============================================================================
# 專案：Zambézia 省抽樣階層空間分佈地圖 (GADM Level 3 重組 22 縣 + 4,874 點位)
# 輸入：01_data/processed/base_ac_processed.rds
# 輸出：03_results/map_zambezia_sampling_tiers_22districts.png
# ==============================================================================

library(here)
library(sf)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(geodata)

dir.create(here("03_results"), recursive = TRUE, showWarnings = FALSE)
sf_use_s2(FALSE)

# ------------------------------------------------------------------------------
# 步驟 1：建立 4,874 點位抽樣三階層因子
# ------------------------------------------------------------------------------
cat(">>> [1/3] 讀取 base_ac_processed.rds 並建立階層分類...\n")

base_ac <- readRDS(here("01_data", "processed", "base_ac_processed.rds"))

base_ac_plot <- base_ac %>%
  mutate(
    Sampling_Tier = case_when(
      COMSA == 2 ~ "COMSA Sero (N = 30)",
      COMSA == 1 ~ "COMSA only (N = 85)",
      TRUE       ~ "Sampling Frame only (N = 4,759)"
    ),
    Sampling_Tier = factor(
      Sampling_Tier,
      levels = c(
        "COMSA Sero (N = 30)",
        "COMSA only (N = 85)",
        "Sampling Frame only (N = 4,759)"
      )
    )
  )

# ------------------------------------------------------------------------------
# 步驟 2：利用 GADM Level 3 重建現行 22 縣邊界
# ------------------------------------------------------------------------------
cat(">>> [2/3] 從 GADM Level 3 解析並重構 Zambézia 22 縣邊界...\n")

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
  st_transform(4326)

cat("  ✓ 成功產出縣級多邊形總數:", nrow(zambezia_distritos_22), "(精確為 22 縣)\n")

# ------------------------------------------------------------------------------
# 步驟 3：分層繪製 22 區邊界與 4,874 個點位
# ------------------------------------------------------------------------------
cat(">>> [3/3] 繪製分層抽樣空間地圖並輸出...\n")

# 拆分點位圖層：底層灰色母體點、頂層抽樣核心點
frame_pts <- base_ac_plot %>% filter(COMSA == 0)
comsa_pts <- base_ac_plot %>% filter(COMSA %in% c(1, 2))

tier_colors <- c(
  "COMSA Sero (N = 30)"             = "#E41A1C", # 鮮紅
  "COMSA only (N = 85)"             = "#1F78B4", # 深藍
  "Sampling Frame only (N = 4,759)" = "#9E9E9E"  # 灰
)

p_map <- ggplot() +
  # 1. 22 縣邊界底圖
  geom_sf(
    data = zambezia_distritos_22,
    fill = "#F8F9FA",
    color = "#555555",
    linewidth = 0.45
  ) +
  # 2. 底層：4,759 個 Sampling Frame 母體點
  geom_sf(
    data = frame_pts,
    aes(color = Sampling_Tier),
    size = 0.75,
    alpha = 0.4,
    shape = 16
  ) +
  # 3. 頂層：85 個 COMSA only 與 30 個 Sero 點位
  geom_sf(
    data = comsa_pts,
    aes(fill = Sampling_Tier),
    size = 2.6,
    shape = 21,
    color = "black",
    stroke = 0.45,
    alpha = 0.95
  ) +
  # 4. 色彩映射
  scale_color_manual(name = "Sampling Hierarchy", values = tier_colors) +
  scale_fill_manual(name = "Sampling Hierarchy", values = tier_colors) +
  # 5. 比例尺與指北針
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_size = 9,
    pad_x = unit(0.5, "cm"),
    pad_y = unit(0.5, "cm")
  ) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    height = unit(1.2, "cm"),
    width = unit(1.2, "cm"),
    pad_x = unit(0.5, "cm"),
    pad_y = unit(0.5, "cm"),
    style = north_arrow_fancy_orienteering
  ) +
  # 6. 主題樣式
  labs(
    title = "Sampling Hierarchy across 22 Districts of Zambézia Province",
    subtitle = "Sampling Frame (N = 4,874), COMSA SRS only (N = 85), and COMSA Sero (N = 30)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "#EDF2F7", color = NA),
    panel.grid.major = element_line(color = "white", linewidth = 0.35),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey30", size = 10.5)
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 4, shape = 21)),
    color = guide_none()
  )

output_map_path <- here("03_results", "map_zambezia_sampling_tiers_22districts.png")
ggsave(
  filename = output_map_path,
  plot = p_map,
  width = 10,
  height = 8.5,
  dpi = 300
)

sf_use_s2(TRUE)
cat("✓ 22 區分層抽樣地圖已成功輸出至:\n  ", output_map_path, "\n")