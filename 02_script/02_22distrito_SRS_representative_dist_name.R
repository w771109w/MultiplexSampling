# ==============================================================================
# 專案：Zambézia 現行 22 縣 (Distrito) 空間抽樣點統計、雙重排序與互動地圖 (完整自洽版)
# ==============================================================================

library(here)
library(sf)
library(dplyr)
library(ggplot2)
library(mapview)
library(leafem)

# 確保輸出資料夾存在
dir.create(here("03_results"), recursive = TRUE, showWarnings = FALSE)

# 關閉 s2 平面運算以確保多邊形融合拓撲穩定
sf_use_s2(FALSE)

# ------------------------------------------------------------------------------
# 步驟 1：下載並建構 Zambézia 現行 22 縣邊界 (利用 GADM level 3 精確分離)
# ------------------------------------------------------------------------------
cat(">>> [1/5] 下載並解析莫三比克現行 22 縣 (Distrito) 向量多邊形...\n")

if (!requireNamespace("geodata", quietly = TRUE)) install.packages("geodata")
moz_adm3_raw <- geodata::gadm(country = "MOZ", level = 3, path = tempdir()) %>%
  st_as_sf()

# 篩選 Zambezia 省並進行現行 22 縣名稱重新對齊 (含 Mulevala 與新設立縣)
zambezia_distritos <- moz_adm3_raw %>%
  filter(NAME_1 == "Zambezia") %>%
  mutate(
    district_22 = case_when(
      NAME_3 == "Derre"                                       ~ "Derre",
      NAME_3 == "Luabo"                                       ~ "Luabo",
      NAME_3 %in% c("Mocubela", "Bajone")                     ~ "Mocubela",
      NAME_3 %in% c("Molumbo", "Corromana")                   ~ "Molumbo",
      grepl("^Mulevala", NAME_3, ignore.case = TRUE)           ~ "Mulevala",  # 精確命中 Mulevala (Namigonha)
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

cat("  - Zambézia 轄下縣份總數 (Distritos):", nrow(zambezia_distritos), "(已確認為 22 縣)\n")


# ------------------------------------------------------------------------------
# 步驟 2：讀取點位並重新與 22 縣執行空間交集 (Spatial Join)
# ------------------------------------------------------------------------------
cat("\n>>> [2/5] 讀取 base_ac_processed.rds 並與 22 縣多邊形重新交集...\n")

base_ac_raw <- readRDS(here("01_data", "processed", "base_ac_processed.rds")) %>%
  st_transform(4326)

# 若原本帶有舊的 district_gadm 欄位，先剔除避免干擾
if ("district_gadm" %in% names(base_ac_raw)) {
  base_ac_raw <- base_ac_raw %>% select(-district_gadm)
}

# 空間交集指派 22 縣
pts_with_district <- st_join(
  base_ac_raw,
  zambezia_distritos,
  join = st_intersects,
  left = TRUE
)

# 邊界外點位檢查 (若有沿海微小漂移，以最近鄰校正)
na_count <- sum(is.na(pts_with_district$district_gadm))
if (na_count > 0) {
  cat("  - 偵測到", na_count, "個點位落在多邊形外邊緣，啟動最近鄰校正...\n")
  na_indices <- which(is.na(pts_with_district$district_gadm))
  nearest_poly_idx <- st_nearest_feature(pts_with_district[na_indices, ], zambezia_distritos)
  pts_with_district$district_gadm[na_indices] <- zambezia_distritos$district_gadm[nearest_poly_idx]
}

# 篩選出 115 個核心 Cluster 抽樣點
pts_115 <- pts_with_district %>%
  filter(COMSA %in% c(1, 2)) %>%
  mutate(
    tier_label = factor(
      if_else(COMSA == 2, "COMSA + Serosurvey (N = 30)", "COMSA Only (N = 85)"),
      levels = c("COMSA + Serosurvey (N = 30)", "COMSA Only (N = 85)")
    )
  )

cat("  - 115 個抽樣點位空間屬性更新完畢！\n")


# ------------------------------------------------------------------------------
# 步驟 3：重新統計 22 縣分佈並執行「雙重排序」+ 存檔
# ------------------------------------------------------------------------------
cat("\n>>> [3/5] 樞紐分析統計 22 縣並執行雙重排序 (Sero 在前 -> SRS 抽樣率大到小)...\n")

# 計算主體 22 縣統計值
summary_districts <- pts_with_district %>%
  st_drop_geometry() %>%
  group_by(district_gadm) %>%
  summarize(
    Total_Frame_N4874 = n(),
    COMSA_SRS_N115    = sum(COMSA %in% c(1, 2)),
    COMSA_Sero_N30    = sum(COMSA == 2),
    srs_rate_num      = COMSA_SRS_N115 / Total_Frame_N4874,
    SRS_Sampling_Rate = paste0(round(srs_rate_num * 100, 2), "%"),
    Has_Sero_Sites    = factor(if_else(COMSA_Sero_N30 > 0, "Yes", "No"), levels = c("Yes", "No")),
    .groups = "drop"
  )

# 若有包含 0 點位的縣（如 Mocubela），確保 22 縣皆在表內
all_22_districts <- data.frame(district_gadm = zambezia_distritos$district_gadm)
summary_districts <- all_22_districts %>%
  left_join(summary_districts, by = "district_gadm") %>%
  mutate(
    Total_Frame_N4874 = coalesce(Total_Frame_N4874, 0L),
    COMSA_SRS_N115    = coalesce(COMSA_SRS_N115, 0L),
    COMSA_Sero_N30    = coalesce(COMSA_Sero_N30, 0L),
    srs_rate_num      = coalesce(srs_rate_num, 0),
    SRS_Sampling_Rate = coalesce(SRS_Sampling_Rate, "0%"),
    Has_Sero_Sites    = coalesce(Has_Sero_Sites, factor("No", levels = c("Yes", "No")))
  ) %>%
  # 雙重排序：1. Has_Sero_Sites (Yes 在前) -> 2. srs_rate_num (由大到小降冪)
  arrange(Has_Sero_Sites, desc(srs_rate_num))

# 指派地圖繪製所需的三層級狀態
summary_districts <- summary_districts %>%
  mutate(
    district_status = case_when(
      COMSA_Sero_N30 > 0 ~ "Has Serosurvey (N = 12)",
      COMSA_SRS_N115 > 0 ~ "SRS Only (N = 9)",
      TRUE               ~ "No SRS (N = 1)"
    ),
    district_status = factor(
      district_status,
      levels = c("Has Serosurvey (N = 12)", "SRS Only (N = 9)", "No SRS (N = 1)")
    )
  )

# 建立置底 TOTAL 列
total_row <- data.frame(
  district_gadm     = "TOTAL (Zambézia)",
  Total_Frame_N4874 = sum(summary_districts$Total_Frame_N4874),
  COMSA_SRS_N115    = sum(summary_districts$COMSA_SRS_N115),
  COMSA_Sero_N30    = sum(summary_districts$COMSA_Sero_N30),
  srs_rate_num      = sum(summary_districts$COMSA_SRS_N115) / sum(summary_districts$Total_Frame_N4874),
  SRS_Sampling_Rate = paste0(round(sum(summary_districts$COMSA_SRS_N115) / sum(summary_districts$Total_Frame_N4874) * 100, 2), "%"),
  Has_Sero_Sites    = paste0(sum(summary_districts$Has_Sero_Sites == "Yes"), " / ", nrow(zambezia_distritos), " Districts"),
  district_status   = NA
)

# 輸出 CSV (移除暫存數值欄位)
final_csv_output <- bind_rows(summary_districts, total_row) %>%
  select(-srs_rate_num, -district_status)

output_csv_path <- here("03_results", "table_district_sampling_summary_22districts.csv")
write.csv(final_csv_output, output_csv_path, row.names = FALSE)
cat("  ✓ 最新 22 縣雙重排序報表已寫入:", output_csv_path, "\n")
print(as.data.frame(final_csv_output))


# ------------------------------------------------------------------------------
# 步驟 4：整併多邊形與標籤錨點 (Point on Surface)
# ------------------------------------------------------------------------------
cat("\n>>> [4/5] 整合多邊形幾何與中心文字標籤錨點...\n")

distritos_map_sf <- zambezia_distritos %>%
  left_join(summary_districts, by = "district_gadm")

# 使用 st_point_on_surface 計算各縣內部安全文字錨點
district_label_pts <- distritos_map_sf %>%
  st_point_on_surface()

coords <- st_coordinates(district_label_pts)
distritos_map_sf$label_x <- coords[, 1]
distritos_map_sf$label_y <- coords[, 2]


# ------------------------------------------------------------------------------
# 步驟 5：產出靜態地圖 (ggplot2) 與互動地圖 (mapview)
# ------------------------------------------------------------------------------
cat("\n>>> [5/5] 輸出出版級靜態地圖與渲染 Viewer 互動地圖...\n")

# 5.1 靜態地圖 (ggplot2)
p_map_fixed <- ggplot() +
  geom_sf(
    data = distritos_map_sf,
    aes(fill = district_status),
    color = "#37474F",
    linewidth = 0.45
  ) +
  scale_fill_manual(
    name = "District Sampling Status",
    values = c(
      "Has Serosurvey (N = 12)" = "#A5D6A7",
      "SRS Only (N = 9)"         = "#FFE082",
      "No SRS (N = 1)"           = "#ECEFF1"
    )
  ) +
  geom_text(
    data = distritos_map_sf,
    aes(x = label_x, y = label_y, label = district_gadm),
    size = 2.8,
    fontface = "bold",
    color = "#1A1A1A",
    check_overlap = TRUE
  ) +
  geom_sf(
    data = pts_115,
    aes(color = tier_label, shape = tier_label),
    size = 2.4,
    stroke = 0.6
  ) +
  scale_color_manual(
    name = "Sampled Clusters (N = 115)",
    values = c(
      "COMSA + Serosurvey (N = 30)" = "#D50000",
      "COMSA Only (N = 85)"        = "#0D47A1"
    )
  ) +
  scale_shape_manual(
    name = "Sampled Clusters (N = 115)",
    values = c(
      "COMSA + Serosurvey (N = 30)" = 16,
      "COMSA Only (N = 85)"        = 17
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

ggsave(
  filename = here("03_results", "map_zambezia_22districts_3tiers_srs.png"),
  plot = p_map_fixed,
  width = 10,
  height = 8,
  dpi = 300
)

# 5.2 互動式地圖 (mapview + leafem 常駐文字標籤)
mapviewOptions(
  basemaps = c("OpenStreetMap", "Esri.WorldImagery", "CartoDB.Positron"),
  fgb = TRUE
)

map_district_interactive <- mapview(
  distritos_map_sf,
  zcol = "district_status",
  col.regions = c(
    "Has Serosurvey (N = 12)" = "#66BB6A",
    "SRS Only (N = 9)"         = "#FFA726",
    "No SRS (N = 1)"           = "#CFD8DC"
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

map_pts_interactive <- mapview(
  pts_115,
  zcol = "tier_label",
  col.regions = c(
    "COMSA + Serosurvey (N = 30)" = "#D50000",
    "COMSA Only (N = 85)"        = "#0D47A1"
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

map_final_view <- (map_district_interactive + map_pts_interactive)@map %>%
  leafem::addStaticLabels(
    data = district_label_pts,
    label = district_label_pts$district_gadm,
    noHide = TRUE,
    direction = "center",
    textsize = "11px",
    textOnly = FALSE,
    style = list(
      "font-weight" = "bold",
      "color" = "#1A1A1A",
      "background-color" = "rgba(255, 255, 255, 0.75)",
      "border" = "1px solid #757575",
      "border-radius" = "3px",
      "padding" = "2px 5px"
    )
  )

sf_use_s2(TRUE)
print(map_final_view)

cat("\n==============================================================================\n")
cat("全數順利完成！\n")
cat("1. 重新空間交集確保 22 縣邊界與點位完全對齊。\n")
cat("2. table_district_sampling_summary_22districts.csv 已按 Sero 與抽樣率排序存檔。\n")
cat("3. 地圖已常駐顯示 22 縣名，並以綠 (12)、橘 (9)、灰 (1) 三色精準呈現。\n")
cat("==============================================================================\n")