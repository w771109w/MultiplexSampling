### calculating tau for COMSA data ###
library(tidyverse)
library(geosphere) 

#-----------Tau Evaluation for single pathogen (bm14)--COMSA data--------------#

# data processing
rm(list=ls())

setwd("/Users/weiliu/Library/CloudStorage/OneDrive-JohnsHopkins/JHU/Project/Moz/01_data")
moz <- read_csv("clean_data_12September2023_STsubset_v2.csv")

###COMSA sample: moz_age

#Rearrange variables.
colnames(moz)
sort(colnames(moz))

moz_analysis <- moz %>% 
  select(sample, age=age_nodecim, gender, district, cluster, area = residence, 
         hhid = hhids_key, lat = latitude, lon = longitude, incom = NationalQuintile, wt = final_wt,
         date = acq51, batch, everything())
colnames(moz_analysis)


#data wrangling
moz_age <- moz_analysis %>% 
  mutate(age_c = case_when(
    age >= 0 & age <=4 ~ "0-59 months",
    age >=5 & age <=17 ~ "5-17 years",
    age >=18 ~"18-49 years",
    age == NA ~ NA
  )) %>%
  mutate(age_c6 = case_when(
    age >= 0 & age <=4 ~ "0-4 years",
    age >=5 & age <=10 ~ "5-10 years",
    age >=11 & age <=20 ~ "11-20 years",
    age >=21 & age <=30 ~ "21-30 years",
    age >=31 & age <=40 ~ "31-40 years",
    age >=41 & age <=49 ~ "41-49 years",
    age == NA ~ NA)) %>%
  mutate(age_sd6 = case_when(
    age >= 0 & age <=4 ~ "0-4 years",
    age >=5 & age <=9 ~ "5-9 years",
    age >=10 & age <=19 ~ "10-19 years",
    age >=20 & age <=29 ~ "20-29 years",
    age >=30 & age <=39 ~ "30-39 years",
    age >=40 & age <=49 ~ "40-49 years",
    age == NA ~ NA)) %>%
  mutate(pos_bm14 = case_when(
    mfi_bm14_107 == "positive" ~ 1,
    mfi_bm14_107 == "negative" ~ 0
  )) %>% 
  mutate(pos_bm33 = case_when(
    mfi_bm33_375 == "positive" ~ 1,
    mfi_bm33_375 == "negative" ~ 0
  )) %>% 
  mutate(pos_wb123 = case_when(
    mfi_wb123_58 == "positive" ~ 1,
    mfi_wb123_58 == "negative" ~ 0
  )) %>% 
  mutate(pos_tmpa = case_when(
    mfi_tmpa_247 == "positive" ~ 1,
    mfi_tmpa_247 == "negative" ~ 0
  )) %>% 
  mutate(pos_rp17 = case_when(
    mfi_rp17_73 == "positive" ~ 1,
    mfi_rp17_73 == "negative" ~ 0
  )) %>% 
  mutate(pos_t24h = case_when(
    mfi_t24h_103 == "positive" ~ 1,
    mfi_t24h_103 == "negative" ~ 0
  )) %>% 
  mutate(pos_es33 = case_when(
    mfi_es33_22 == "positive" ~ 1,
    mfi_es33_22 == "negative" ~ 0
  )) %>% 
  mutate(pos_sm25 = case_when(
    mfi_sm25_38 == "positive" ~ 1,
    mfi_sm25_38 == "negative" ~ 0
  )) %>% 
  mutate(pos_sea = case_when(
    mfi_sea_410 == "positive" ~ 1,
    mfi_sea_410 == "negative" ~ 0
  )) %>% 
  mutate(pos_pgp3 = case_when(
    mfi_pgp3_212 == "positive" ~ 1,
    mfi_pgp3_212 == "negative" ~ 0
  )) %>% 
  mutate(pos_ct694 = case_when(
    mfi_ct694_108 == "positive" ~ 1,
    mfi_ct694_108 == "negative" ~ 0
  )) %>% 
  mutate(pos_nie = case_when(
    mfi_nie_526 == "positive" ~ 1,
    mfi_nie_526 == "negative" ~ 0
  )) %>% 
  mutate(pos_ov16 = case_when(
    mfi_ov16_422 == "positive" ~ 1,
    mfi_ov16_422 == "negative" ~ 0
  )) %>% 
  mutate(pos_vsp5 = case_when(
    mfi_vsp5_148 == "positive" ~ 1,
    mfi_vsp5_148 == "negative" ~ 0
  )) %>% 
  mutate(pos_vsp3 = case_when(
    mfi_vsp3_72 == "positive" ~ 1,
    mfi_vsp3_72 == "negative" ~ 0
  )) %>% 
  mutate(pos_cp23 = case_when(
    mfi_cp23_830 == "positive" ~ 1,
    mfi_cp23_830 == "negative" ~ 0
  )) %>% 
  mutate(pos_cp17 = case_when(
    mfi_cp17_72 == "positive" ~ 1,
    mfi_cp17_72 == "negative" ~ 0
  )) %>% 
  mutate(mfi_bm14x = case_when(mfi_bm14 <= 0 ~ 1, #assign 1 to 0, ensure no NaN after log
                               TRUE ~ as.numeric(mfi_bm14))) %>%
  mutate(log_bm14 = log10(mfi_bm14x)) %>%
  mutate(mfi_bm33x = case_when(mfi_bm33 <= 0 ~ 1,
                               TRUE ~ as.numeric(mfi_bm33))) %>%
  mutate(log_bm33 = log10(mfi_bm33x)) %>%
  mutate(mfi_wb123x = case_when(mfi_wb123 <= 0 ~ 1,
                                TRUE ~ as.numeric(mfi_wb123))) %>%
  mutate(log_wb123 = log10(mfi_wb123x)) %>%
  mutate(mfi_tmpax = case_when(mfi_tmpa <= 0 ~ 1,
                               TRUE ~ as.numeric(mfi_tmpa))) %>%
  mutate(log_tmpa = log10(mfi_tmpax)) %>%
  mutate(mfi_rp17x = case_when(mfi_rp17 <= 0 ~ 1,
                               TRUE ~ as.numeric(mfi_rp17))) %>%
  mutate(log_rp17 = log10(mfi_rp17x)) %>%
  mutate(mfi_t24hx = case_when(mfi_t24h <= 0 ~ 1,
                               TRUE ~ as.numeric(mfi_t24h))) %>%
  mutate(log_t24h = log10(mfi_t24hx)) %>%
  mutate(mfi_es33x = case_when(mfi_es33 <= 0 ~ 1,
                               TRUE ~ as.numeric(mfi_es33))) %>%
  mutate(log_es33 = log10(mfi_es33x)) %>%
  mutate(mfi_sm25x = case_when(mfi_sm25 <= 0 ~ 1,
                               TRUE ~ as.numeric(mfi_sm25))) %>%
  mutate(log_sm25 = log10(mfi_sm25x)) %>%
  mutate(mfi_seax = case_when(mfi_sea <= 0 ~ 1,
                              TRUE ~ as.numeric(mfi_sea))) %>%
  mutate(log_sea = log10(mfi_seax)) %>%
  mutate(mfi_pgp3x = case_when(mfi_pgp3 <= 0 ~ 1,
                               TRUE ~ as.numeric(mfi_pgp3))) %>%
  mutate(log_pgp3 = log10(mfi_pgp3x)) %>%
  mutate(mfi_ct694x = case_when(mfi_ct694 <= 0 ~ 1,
                                TRUE ~ as.numeric(mfi_ct694))) %>%
  mutate(log_ct694 = log10(mfi_ct694x)) %>%
  mutate(mfi_niex = case_when(mfi_nie <= 0 ~ 1,
                              TRUE ~ as.numeric(mfi_nie))) %>%
  mutate(log_nie = log10(mfi_niex)) %>%
  mutate(mfi_ov16x = case_when(mfi_ov16 <= 0 ~ 1,
                               TRUE ~ as.numeric(mfi_ov16))) %>%
  mutate(log_ov16 = log10(mfi_ov16x)) %>%
  mutate(mfi_vsp5x = case_when(mfi_vsp5 <= 0 ~ 1,
                               TRUE ~ as.numeric(mfi_vsp5))) %>%
  mutate(log_vsp5 = log10(mfi_vsp5x)) %>%
  mutate(mfi_vsp3x = case_when(mfi_vsp3 <= 0 ~ 1,
                               TRUE ~ as.numeric(mfi_vsp3))) %>%
  mutate(log_vsp3 = log10(mfi_vsp3x)) %>%
  mutate(mfi_cp23x = case_when(mfi_cp23 <= 0 ~ 1,
                               TRUE ~ as.numeric(mfi_cp23))) %>%
  mutate(log_cp23 = log10(mfi_cp23x)) %>%
  mutate(mfi_cp17x = case_when(mfi_cp17 <= 0 ~ 1,
                               TRUE ~ as.numeric(mfi_cp17))) %>%
  mutate(log_cp17 = log10(mfi_cp17x)) %>%
  mutate(overall = 1) %>%
  mutate(mfi_gst_cdc_x = case_when(mfi_gst_cdc <= 0 ~ 1,
                                   TRUE ~ as.numeric(mfi_gst_cdc)))%>%
  mutate(mfi_gst_ls_x = case_when(mfi_gst_ls <= 0 ~ 1,
                                  TRUE ~ as.numeric(mfi_gst_ls)))%>%
  mutate(log_gst_cdc = log(mfi_gst_cdc_x),
         log_gst_ls = log(mfi_gst_ls_x))

rm(moz_analysis)
rm(moz)




# build cluster's polygon and convert it to sf object





# aggregate data by cluster
cluster_data <- moz_age %>%
  group_by(cluster) %>%
  summarise(
    s = sum(overall),
    w = sum(wt),
    p = mean(pos_bm14),
    lon = mean(lon),
    lat = mean(lat)
  )


### tau functions for survey data ###

# install.packages("geosphere") # for calculating distance matrix in km by lon & lat

# check distance matrix (unit: km)
distance_matrix <- as.matrix(distm(cluster_data[, c("lon", "lat")], fun = distHaversine))/1000
max_value <- max(distance_matrix)


distance_matrix <- as.matrix(distm(moz_age[, c("lon", "lat")], fun = distHaversine))/1000
max_value <- max(distance_matrix)
summary(distance_matrix[1,])
head(distance_matrix[1,],50)



####------------ tau function (excluding self) ------------#### 
calculate_tau_r_excl <- function(data, r2) {
  library(geosphere) 
  
  # calculate distance matrix (unit: km)
  distance_matrix <- as.matrix(distm(data[, c("lon", "lat")], fun = distHaversine))/1000
  
  # population weighted seroprevalence Pi(.)
  total_weighted_infections <- sum(data$s * data$w * data$p)
  total_weighted_samples <- sum(data$s * data$w)
  Pi_all <- total_weighted_infections / total_weighted_samples
  
  # create empty vector for pi_i(r2)
  pi_values <- numeric(nrow(data))
  
  # compute pi_i(r2) for all clusters
  for (i in 1:nrow(data)) {
    # find cluster within (0, r2] (exclude itself)
    neighbors <- which(distance_matrix[i, ] > 0 & distance_matrix[i, ] <= r2)
    
    # skip if no neighbors
    if (length(neighbors) == 0) next
    
    # compute pi_i(r2)
    s_w_p_j <- sum(data$s[neighbors] * data$w[neighbors] * data$p[neighbors])
    s_w_j <- sum(data$s[neighbors] * data$w[neighbors])
    pi_values[i] <- s_w_p_j / s_w_j
  }
  
  # compute Pi(r2)
  valid_indices <- which(pi_values > 0)
  Pi_r <- sum(pi_values[valid_indices] * data$s[valid_indices] * data$w[valid_indices] * data$p[valid_indices]) /
    total_weighted_infections
  
  # compute tau(r2)
  tau_value <- Pi_r / Pi_all
  
  return(list(Pi_r = Pi_r, Pi_all = Pi_all, tau = tau_value))
}


####------------ tau function (including self) ------------#### 
calculate_tau_r_incl <- function(data, r2) {
  library(geosphere)
  
  # calculate distance matrix (unit: km)
  distance_matrix <- as.matrix(distm(data[, c("lon", "lat")], fun = distHaversine))/1000
  
  # population weighted seroprevalence Pi(.)
  total_weighted_infections <- sum(data$s * data$w * data$p)
  total_weighted_samples <- sum(data$s * data$w)
  Pi_all <- total_weighted_infections / total_weighted_samples
  
  # create empty vector for pi_i(r2)
  pi_values <- numeric(nrow(data))
  
  # compute pi_i(r2) for all clusters
  for (i in 1:nrow(data)) {
    # find cluster within (0, r2] (include itself)
    neighbors <- which(distance_matrix[i, ] >= 0 & distance_matrix[i, ] <= r2)
    
    # skip if no neighbors
    if (length(neighbors) == 0) next
    
    # compute pi_i(r2)
    s_w_p_j <- sum(data$s[neighbors] * data$w[neighbors] * data$p[neighbors])
    s_w_j <- sum(data$s[neighbors] * data$w[neighbors])
    pi_values[i] <- s_w_p_j / s_w_j
  }
  
  # compute Pi(r2)
  valid_indices <- which(pi_values > 0)
  Pi_r <- sum(pi_values[valid_indices] * data$s[valid_indices] * data$w[valid_indices] * data$p[valid_indices]) /
    total_weighted_infections
  
  # compute tau(r2)
  tau_value <- Pi_r / Pi_all
  
  return(list(Pi_r = Pi_r, Pi_all = Pi_all, tau = tau_value))
}




####------------  tau function (excluding self, without weight) ------------#### 
calculate_tau_r_excl_no_wt <- function(data, r2) {
  library(geosphere) 
  
  # calculate distance matrix (unit: km)
  distance_matrix <- as.matrix(distm(data[, c("lon", "lat")], fun = distHaversine))/1000
  
  # population weighted seroprevalence Pi(.)
  total_infections <- sum(data$s * data$p)
  total_samples <- sum(data$s)
  Pi_all <- total_infections / total_samples
  
  # create empty vector for pi_i(r2)
  pi_values <- numeric(nrow(data))
  
  # compute pi_i(r2) for all clusters
  for (i in 1:nrow(data)) {
    # find cluster within (0, r2] (exclude itself)
    neighbors <- which(distance_matrix[i, ] > 0 & distance_matrix[i, ] <= r2)
    
    # skip if no neighbors
    if (length(neighbors) == 0) next
    
    # compute pi_i(r2)
    s_p_j <- sum(data$s[neighbors] * data$p[neighbors])
    s_j <- sum(data$s[neighbors])
    pi_values[i] <- s_p_j / s_j
  }
  
  # compute Pi(r2)
  valid_indices <- which(pi_values > 0)
  Pi_r <- sum(pi_values[valid_indices] * data$s[valid_indices] * data$p[valid_indices]) /
    total_infections
  
  # compute tau(r2)
  tau_value <- Pi_r / Pi_all
  
  return(list(Pi_r = Pi_r, Pi_all = Pi_all, tau = tau_value))
}



####------------   tau function (including self, without weight) ------------####
calculate_tau_r_incl_no_wt <- function(data, r2) {
  library(geosphere) 
  
  # calculate distance matrix (unit: km)
  distance_matrix <- as.matrix(distm(data[, c("lon", "lat")], fun = distHaversine))/1000
  
  # population weighted seroprevalence Pi(.)
  total_infections <- sum(data$s * data$p)
  total_samples <- sum(data$s)
  Pi_all <- total_infections / total_samples
  
  # create empty vector for pi_i(r2)
  pi_values <- numeric(nrow(data))
  
  # compute pi_i(r2) for all clusters
  for (i in 1:nrow(data)) {
    # find cluster within (0, r2] (include itself)
    neighbors <- which(distance_matrix[i, ] >= 0 & distance_matrix[i, ] <= r2)
    
    # skip if no neighbors
    if (length(neighbors) == 0) next
    
    # compute pi_i(r2)
    s_p_j <- sum(data$s[neighbors] * data$p[neighbors])
    s_j <- sum(data$s[neighbors])
    pi_values[i] <- s_p_j / s_j
  }
  
  # compute Pi(r2)
  valid_indices <- which(pi_values > 0)
  Pi_r <- sum(pi_values[valid_indices] * data$s[valid_indices] * data$p[valid_indices]) /
    total_infections
  
  # compute tau(r2)
  tau_value <- Pi_r / Pi_all
  
  return(list(Pi_r = Pi_r, Pi_all = Pi_all, tau = tau_value))
}



####------------ Running for bm14 only as an example ------------####

# calculate tau for different r values
r_vector <- seq(0, max(distance_matrix)/3, by = 10)
tau_results_excl <- sapply(r_vector, function(r) calculate_tau_r_excl(cluster_data, r)$tau)
tau_results_incl <- sapply(r_vector, function(r) calculate_tau_r_incl(cluster_data, r)$tau)
tau_results_excl_no_wt <- sapply(r_vector, function(r) calculate_tau_r_excl_no_wt(cluster_data, r)$tau)
tau_results_incl_no_wt <- sapply(r_vector, function(r) calculate_tau_r_incl_no_wt(cluster_data, r)$tau)

# plot tau vs distance
# dev.off()  # Closes any open graphics devices
tau_data <- 
  data.frame(distance = r_vector, tau_excl = tau_results_excl,
             tau_incl = tau_results_incl, tau_excl_no_wt = tau_results_excl_no_wt, 
             tau_incl_no_wt = tau_results_incl_no_wt)

library(ggplot2)
ggplot(data = tau_data, aes(x= distance)) +
  geom_line(aes(y = tau_excl, color = "Excluding Self")) +
  geom_point(aes(y = tau_excl, color = "Excluding Self")) +
  geom_line(aes(y = tau_incl, color = "Including Self")) +
  geom_point(aes(y = tau_incl, color = "Including Self")) +
  geom_line(aes(y = tau_excl_no_wt, color = "Excluding Self (No Wt)")) +
  geom_point(aes(y = tau_excl_no_wt, color = "Excluding Self (No Wt)")) +
  geom_line(aes(y = tau_incl_no_wt, color = "Including Self (No Wt)")) +
  geom_point(aes(y = tau_incl_no_wt, color = "Including Self (No Wt)")) +
  geom_hline( yintercept = 1,linetype = "dashed") +
  scale_color_manual(
    name = "Tau Method",values = c("Excluding Self" = "#1ABC9C", 
                                   "Excluding Self (No Wt)" = "#154360",
                                   "Including Self" = "#F39B7FB2",
                                   "Including Self (No Wt)" = "#DC0000B2")) +
  labs(x = "Distance (km)", y = "Tau", title = "Tau(bm14) vs Distance") +
  theme_minimal() +
  theme(legend.position = c(0.8, 0.9))

####------------ Running for all antigens ------------####


moz_selected <- moz_age %>% 
  select(matches("pos_"), wt, cluster, lon, lat)
pos_list <- colnames(moz_selected |> select(matches("pos_")))


cluster_list <- lapply(pos_list, function(pos_col){
  moz_selected <- moz_age |> select(matches(pos_col), wt, cluster, lon, lat)
  
  # aggregate data by cluster
  cluster_data <- moz_selected %>%
    group_by(cluster) %>%
    summarise(
      s = n(),
      w = sum(wt),
      p = mean(!!sym(pos_col)),
      lon = mean(lon),
      lat = mean(lat)
    )
  
  # calculate tau for different r values
  r_vector <- seq(0, max(distance_matrix)/3, by = 10)
  tau_results_excl <- sapply(r_vector, function(r) calculate_tau_r_excl(cluster_data, r)$tau)
  tau_results_incl <- sapply(r_vector, function(r) calculate_tau_r_incl(cluster_data, r)$tau)
  tau_results_excl_no_wt <- sapply(r_vector, function(r) calculate_tau_r_excl_no_wt(cluster_data, r)$tau)
  tau_results_incl_no_wt <- sapply(r_vector, function(r) calculate_tau_r_incl_no_wt(cluster_data, r)$tau)
  
  # plot tau vs distance
  # dev.off()  # Closes any open graphics devices
  tau_data <- 
    data.frame(distance = r_vector, tau_excl = tau_results_excl,
               tau_incl = tau_results_incl, tau_excl_no_wt = tau_results_excl_no_wt, 
               tau_incl_no_wt = tau_results_incl_no_wt)
  
  return(tau_data)
  
})

####------------ Plotting for all antigens ------------####

library(ggplot2)
library(purrr) 

names(cluster_list) <- pos_list

plots <- map2(cluster_list, names(cluster_list), function(data, name) {
  ggplot(data, aes(x = distance)) +
    geom_line(aes(y = tau_incl, color = "Including Self"), size = 1.2) +  # Thicker lines
    geom_point(aes(y = tau_incl, color = "Including Self"), size = 2) +  # Larger points
    geom_line(aes(y = tau_incl_no_wt, color = "Including Self (No Wt)"), size = 1.2) +
    geom_point(aes(y = tau_incl_no_wt, color = "Including Self (No Wt)"), size = 2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_color_manual(
      name = "Tau Method",
      values = c(
        "Including Self" = "#F39B7FB2",
        "Including Self (No Wt)" = "#DC0000B2"
      )
    ) +
    scale_y_continuous(limits = c(0, 3.5), breaks = seq(0, 3.5, by = 0.5)) +  # Set y-axis intervals to 1
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +  # Set x-axis intervals to 10
    labs(x = "Distance (km)", y = "Tau", title = paste("Tau(", name, ")")) +
    theme_classic() +
    theme(legend.position = "none")
})

print(plots)

library(patchwork)

combined_plot <- wrap_plots(plots, ncol = 3)

print(combined_plot)
