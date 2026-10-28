library(sf)
library(tidyverse)
library(ggsflabel)
library(ggspatial)
library(ggOceanMaps)
library(patchwork)

setwd("~/Dropbox/COMSA - Serosurveillance/Data/Saki/")

## Read in COMSA cluster spatial data
dat_COMSA <- st_read("comsa.zambezia_with_clusters/comsa.zambezia_with_clusters.shp")

## Pull in shapefile for Mozambique
## From: https://data.humdata.org/dataset/cod-ab-moz
Mozambique <- st_read("moz_adm_20190607b_shp/moz_admbnda_adm1_ine_20190607.shp")

## Extract Zambezia province only
Mozambique %>%
  filter(ADM1_PT == "Zambezia") -> Zambezia

## Plot just the 30 DBS clusters
dat_COMSA %>%
  filter(sero==1) %>%
  group_by(cluster) %>%
  summarize(geometry = st_union(geometry)) %>% 
  ungroup() %>%
  ggplot() +
  geom_sf(data=Zambezia, size=0.05, fill="lightgreen", colour="black") +
  geom_sf(size=0.1, fill="black") +
  ## Add cluster label
  geom_sf_label_repel(data=dat_COMSA %>% filter(sero==1) %>% group_by(cluster) %>% slice(1) %>% ungroup(), aes(label=cluster), force=30, size=3, label.size=0.1, segment.size=0.1, colour="red") +
  theme_classic() +
  #ggtitle("30 DBS clusters in Zambezia province, Mozambique") +
  # annotation_scale(location = "br", 
  #                  width_hint = 0.25, 
  #                  text_cex = 1) +
  # annotation_north_arrow(location = "br", 
  #                        which_north = "true", 
  #                        pad_x = unit(0.15, "in"),
  #                        pad_y = unit(0.3, "in"),
  #                        style = north_arrow_fancy_orienteering) +
  xlab("") + ylab("") -> plot_map_Zambezia

## Fancier map, zoomed out
basemap(data=Mozambique, land.col="white", bathymetry=FALSE) +
  geom_sf(data=Mozambique) +
  geom_sf(data=Zambezia, size=4, fill="lightgreen") +
  geom_sf(data=dat_COMSA %>% filter(sero==1), fill="black") +
  theme_minimal() +
  ggtitle("Zambezia province, Mozambique") +
  annotation_scale(location = "br", 
                   width_hint = 0.25, 
                   text_cex = 1) +
  annotation_north_arrow(location = "br", 
                         which_north = "true", 
                         pad_x = unit(0.15, "in"),
                         pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  xlab("") + ylab("") +
  theme(panel.background = element_rect(fill = "lightblue"),
        panel.ontop = FALSE,
        panel.grid = element_blank()) -> plot_map_Moz

plot_map_Moz + plot_map_Zambezia
