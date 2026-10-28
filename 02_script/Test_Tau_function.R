
##### Tau function for survey data #####



clusters <- data.frame(
  id = 1:5,
  x = c(1, 2, 5, 8, 9),  # X 
  y = c(1, 3, 5, 8, 10), # Y 
  s = c(100, 200, 150, 300, 250), # sample size
  w = c(1.2, 0.8, 1.5, 0.9, 1.1), # weight
  p = c(0.10, 0.15, 0.20, 0.12, 0.18) # Seroprevalence
)
clusters <- clusters |> arrange(s)

data <- clusters

distance_matrix <- as.matrix(dist(clusters[, c("x", "y")]))
pi_values <- numeric(nrow(clusters))

total_weighted_infections <- sum(data$s * data$w * data$p)
total_weighted_samples <- sum(data$s * data$w)
Pi_all <- total_weighted_infections / total_weighted_samples


library(dplyr)


calculate_tau_r_incl <- function(data, r2) {

  # 計算距離矩陣
  distance_matrix <- as.matrix(dist(data[, c("x", "y")]))
  
  # 計算全體加權 seroprevalence Pi(.)
  total_weighted_infections <- sum(data$s * data$w * data$p)
  total_weighted_samples <- sum(data$s * data$w)
  Pi_all <- total_weighted_infections / total_weighted_samples
  
  # 初始化存儲 pi_i(r2) 的向量
  pi_values <- numeric(nrow(data))
  
  # 遍歷所有 clusters，計算 pi_i(r2)
  for (i in 1:nrow(data)) {
    # 找出距離在 (0, r2] 內的 clusters (包含自己)
    neighbors <- unname(which(distance_matrix[i, ] >= 0 & distance_matrix[i, ] <= r2))
    
    # 若沒有鄰居，則跳過
    if (length(neighbors) == 0) next
    
    # 計算 pi_i(r2)
    s_w_p_j <- sum(data$s[neighbors] * data$w[neighbors] * data$p[neighbors])
    s_w_j <- sum(data$s[neighbors] * data$w[neighbors])
    pi_values[i] <- s_w_p_j / s_w_j
  }
  
  # 計算 Pi(r2)
  valid_indices <- which(pi_values > 0) # 只考慮有數據的 clusters
  Pi_r <- sum(pi_values[valid_indices] * data$s[valid_indices] * data$w[valid_indices] * data$p[valid_indices]) /
    total_weighted_infections
  
  # 計算 tau(r2)
  tau_value <- Pi_r / Pi_all
  
  return(list(Pi_r = Pi_r, Pi_all = Pi_all, tau = tau_value))
}




calculate_tau_r_excl <- function(data, r2) {
  
  # 計算距離矩陣
  distance_matrix <- as.matrix(dist(data[, c("x", "y")]))
  
  # 計算全體加權 seroprevalence Pi(.)
  total_weighted_infections <- sum(data$s * data$w * data$p)
  total_weighted_samples <- sum(data$s * data$w)
  Pi_all <- total_weighted_infections / total_weighted_samples
  
  # 初始化存儲 pi_i(r2) 的向量
  pi_values <- numeric(nrow(data))
  
  # 遍歷所有 clusters，計算 pi_i(r2)
  for (i in 1:nrow(data)) {
    # 找出距離在 (0, r2] 內的 clusters (不包含自己)
    neighbors <- unname(which(distance_matrix[i, ] > 0 & distance_matrix[i, ] <= r2))
    
    # 若沒有鄰居，則跳過
    if (length(neighbors) == 0) next
    
    # 計算 pi_i(r2)
    s_w_p_j <- sum(data$s[neighbors] * data$w[neighbors] * data$p[neighbors])
    s_w_j <- sum(data$s[neighbors] * data$w[neighbors])
    pi_values[i] <- s_w_p_j / s_w_j
  }
  
  # 計算 Pi(r2)
  valid_indices <- which(pi_values > 0) # 只考慮有數據的 clusters
  Pi_r <- sum(pi_values[valid_indices] * data$s[valid_indices] * data$w[valid_indices] * data$p[valid_indices]) /
    total_weighted_infections
  
  # 計算 tau(r2)
  tau_value <- Pi_r / Pi_all
  
  return(list(Pi_r = Pi_r, Pi_all = Pi_all, tau = tau_value))
}



r_vector <- seq(min(distance_matrix[,1]), max(distance_matrix[,1])+1, by = 1)
tau_results_incl <- sapply(r_vector, function(r) calculate_tau_r_incl(clusters, r)$tau)
tau_results_excl <- sapply(r_vector, function(r) calculate_tau_r_excl(clusters, r)$tau)


ggplot() +
  geom_line(aes(x = r_vector, y = tau_results_incl), color = "blue", size = 1) +
  geom_line(aes(x = r_vector, y = tau_results_excl), color = "red", size = 1) +
  labs(title = "Tau function", x = "r", y = "Tau(r)") +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
  theme_classic()
