####Bootstrap 
library(MASS)
source('R_code/Utility_func.R')
set.seed(10212024)

###_______Functions__________###
average_distance <- function(distance_matrix, idx_celltype,Frequency,D_Unknown) {
  # Extract specific nodes based on the provided sample indices
  cell_type <- rownames(distance_matrix)[idx_celltype]
  Freq<-Frequency
  repetition_df <- data.frame(cell_type,Freq)
  
  node_count <- aggregate(Freq ~ cell_type, data = repetition_df, sum)
  unique_node <- node_count$cell_type
  
  # Initialize variable to store the average distance
  N=20
  weighted_average_distance <- NA
  
  # Handle cases with only one unique node
  if (length(unique_node) == 1) {
    if (unique_node == 'Unknown') {
      weighted_average_distance <- D_Unknown  # Special case for 'Unknown'
    } else {
      weighted_average_distance <- 0  # If only one node and not 'Unknown'
    }
  } else {
    # Subset the distance matrix for the specific nodes, excluding "Unknown" if it exists
    filtered_nodes <- unique_node[unique_node != "Unknown"]
    
    if (length(filtered_nodes) > 0) {
      specific_distance_matrix <- distance_matrix[filtered_nodes, filtered_nodes, drop = FALSE]
      
      # Create a named vector for repetition counts
      specific_repetition_counts <- setNames(repetition_df$Freq, repetition_df$cell_type)
      
      # Initialize total weighted distance and total weight
      total_weighted_distance <- 0
      total_weight <- 0.5 * (N - 1) * N  # Based on the formula
      
      # Calculate the weighted distance for known nodes
      for (i in 1:nrow(specific_distance_matrix)) {
        for (j in i:ncol(specific_distance_matrix)) {
          # Get the pairwise distance between nodes i and j
          distance <- specific_distance_matrix[i, j]
          
          # Get the repetition counts for nodes i and j using node names
          rep_i <- specific_repetition_counts[filtered_nodes[i]]
          rep_j <- specific_repetition_counts[filtered_nodes[j]]
          
          # Calculate the weighted distance
          weighted_distance <- distance * rep_i * rep_j
          
          # Add the weighted distance to the total
          total_weighted_distance <- total_weighted_distance + weighted_distance
        }
      }
      
      # Calculate the weighted distance for "Unknown" nodes if they exist
      n_unknown <- specific_repetition_counts["Unknown"]
      D_unknown <- D_Unknown
      if (!is.na(n_unknown) && n_unknown > 0) {
        # Add the distance for "Unknown" nodes (self-distances)
        total_distance <- total_weighted_distance + (D_unknown * 0.5 * n_unknown * (n_unknown - 1))
        
        # Add distances between "Unknown" nodes and all other nodes
        for (i in 1:length(filtered_nodes)) {
          rep_i <- specific_repetition_counts[filtered_nodes[i]]
          total_distance <- total_distance + (D_unknown * rep_i * n_unknown)
        }
      } else {
        total_distance <- total_weighted_distance
      }
      
      # Calculate the weighted average distance
      weighted_average_distance <- as.numeric(total_distance / total_weight)
    }
  }
  
  return(weighted_average_distance)  # Return the average distance
}


bootstrap_average_distance <- function(distance_matrix, n_run, sample_size,Frequency,D_Unknown) {
  # Initialize a vector to store the results of each bootstrap iteration
  bootstrap_results <- numeric(n_run)
  
  # Loop for each bootstrap sample
  for (i in 1:n_run) {
    # Randomly sample 5 cell types with replacement from the 35 available
    sampled_indices <- sample(1:nrow(distance_matrix), sample_size, replace = TRUE)
    
    # Calculate the average distance using the average_distance function
    bootstrap_results[i] <- average_distance(distance_matrix, sampled_indices,Frequency,D_Unknown)
  }
  
  return(bootstrap_results)
}


###________Bootstrap____###

bootstrap_res1_10k<-bootstrap_average_distance(pbmc_dis_mtx_na,n_run = 10000, sample_size = 4,Frequency=c(1,8,1,10),D_Unknown=117.3096)
bootstrap_res2_10k<-bootstrap_average_distance(pbmc_dis_mtx_na,n_run = 10000, sample_size = 5,Frequency=c(2,12,1,1,4),D_Unknown=117.3096)
bootstrap_res3_10k<-bootstrap_average_distance(pbmc_dis_mtx_na,n_run = 10000, sample_size = 6,Frequency=c(2,11,1,2,1,3),D_Unknown=117.3096)
bootstrap_res4_10k<-bootstrap_average_distance(pbmc_dis_mtx_na,n_run = 10000, sample_size = 5,Frequency=c(1,5,1,9,4),D_Unknown=117.3096)
bootstrap_res5_10k<-bootstrap_average_distance(pbmc_dis_mtx_na,n_run = 10000, sample_size = 4,Frequency=c(6,5,3,6),D_Unknown=117.3096)
bootstrap_res6_10k<-bootstrap_average_distance(pbmc_dis_mtx_na,n_run = 10000, sample_size = 5,Frequency=c(5,4,3,2,6),D_Unknown=117.3096)
bootstrap_res7_10k<-bootstrap_average_distance(pbmc_dis_mtx_na,n_run = 10000, sample_size = 4,Frequency=c(1,4,1,14),D_Unknown=117.3096)
bootstrap_res8_10k<-bootstrap_average_distance(pbmc_dis_mtx_na,n_run = 10000, sample_size = 3 ,Frequency=c(5,8,7),D_Unknown=117.3096)
bootstrap_res9_10k<-bootstrap_average_distance(pbmc_dis_mtx_na,n_run = 10000, sample_size = 3 ,Frequency=c(18,2,1),D_Unknown=117.3096)

boot_df<-data.frame(cbind(bootstrap_res1_10k,bootstrap_res2_10k,bootstrap_res3_10k,
                          bootstrap_res4_10k,bootstrap_res5_10k,bootstrap_res6_10k,
                          bootstrap_res7_10k,bootstrap_res8_10k,bootstrap_res9_10k
                          ))
#
#hist(bootstrap_res1_10k,breaks=20,xlim =c(0,80),main=c('Freq=c(1,8,1,10)'))
#hist(bootstrap_res2_10k,breaks=20,xlim =c(0,120),main=c('Freq=c(2,12,1,1,4)'))


#hist(log(bootstrap_res_10k),breaks=20,xlim = c(-4,6),main=c('Freq=c(1,8,1,10)'))
#hist(log(bootstrap_res2_10k),breaks=20,xlim =c(-4,6),main=c('Freq=c(2,12,1,1,4)'))

###Fit function
fit_plot<-function(fit_data,freq){
  par(mfrow=c(1,2))
  ###log normal fit
  fit_data<-fit_data[fit_data > 0]
  normal_fit <- fitdistr(log(fit_data), "normal")
  fitted_mean <- normal_fit$estimate["mean"]
  fitted_sd <- normal_fit$estimate["sd"]
  
  hist(log(fit_data), breaks = 40, probability = TRUE, main = c(paste0('freq=',freq,'normal fit')),
       xlab = "log(data)", ylab = "Density", col = "gray", border = "black",xlim = c(-4,5))
  curve(dnorm(x, mean = fitted_mean, sd = fitted_sd), col = "blue", lwd = 2, add = TRUE)
  
  ### Gamma fit
  gamma_fit <- glm(boot_all_nozero ~ 1, family = Gamma(link = 'log'))
  fitted_mean_g <- exp(coef(gamma_fit))
  gamma_shape <- 1 / summary(gamma_fit)$dispersion
  gamma_rate <- gamma_shape / fitted_mean_g
  
  hist(fit_data, breaks = 40, probability = TRUE, xlim = c(0, 120), 
       main = paste0('freq=', freq, ' Gamma fit'), ylim = c(0, 0.06),
       xlab = "Data", ylab = "Density", col = "gray", border = "black")
  curve(dgamma(x, shape = gamma_shape, rate = gamma_rate), col = "blue", lwd = 2, add = TRUE)
  
}  


###negtive binomal
  nb_fit <- glm.nb(fit_data ~ 1)
  fitted_mean_NB <- exp(coef(nb_fit))
  theta <- nb_fit$theta
  x_vals <- 0:max(fit_data)
  fitted_density <- dnbinom(x_vals, size = theta, mu = fitted_mean_NB)
  hist(fit_data,breaks=30,xlim =c(0,120),main=c(paste0('freq=',freq,'NB fit')),probability = TRUE,ylim=c(0,0.06))
  lines(x_vals, fitted_density, col = "blue", lwd = 2)
  
###Inverse gamma
  inverse_data<- 1 / fit_data
  inverse_gamma_fit <- glm( inverse_data ~ 1, family = Gamma)
  fitted_mean_ig <-  exp(coef(inverse_gamma_fit))
  ig_shape <- 1 / summary(inverse_gamma_fit)$dispersion
  ig_rate <- ig_shape / fitted_mean_ig
  
  hist(fit_data, breaks = 30, probability = TRUE, xlim = c(0, 80), 
       main = paste0('freq=', freq, ' Inverse Gamma fit'), ylim = c(0, 0.06),
       xlab = "Data", ylab = "Density", col = "gray", border = "black")
  curve(dgamma(x, shape = ig_shape, rate = ig_rate) , col = "blue", lwd = 2, add = TRUE)


####Combine the 9 scenarios
boot_all<-c(bootstrap_res1_10k,bootstrap_res2_10k,bootstrap_res3_10k,
            bootstrap_res4_10k,bootstrap_res5_10k,bootstrap_res6_10k,
            bootstrap_res7_10k,bootstrap_res8_10k,bootstrap_res9_10k)
fit_plot(boot_all,freq='combine')


###Gamma fit using fitdistr function
fitdistr(x, "gamma", start=list(shape=1, rate=1))$estimate

boot_all_nozero<-boot_all[boot_all>0]
gamma_fit <- fitdistr(boot_all_nozero, "gamma",start=list(shape=1, rate=1))

# Extract the estimated parameters
shape <- gamma_fit $estimate["shape"]
rate <- gamma_fit $estimate["rate"]


###Histogram showing the distribution
par(mfrow=c(1,1))
hist(boot_all, breaks = 40, probability = TRUE, xlim = c(0, 80), 
     main =  'Gamma fit', ylim = c(0, 0.06),
     xlab = "Data", ylab = "Density", col = "gray", border = "white")
curve(dgamma(x, shape = shape, rate = rate), col = "blue", lwd = 2, add = TRUE)

# Sample data
x <- 1:10
y <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29)

# Basic line plot
plot(dgamma(x, shape = shape, rate = rate), type = "o", col = "blue", xlab = "X-axis Label", ylab = "Y-axis Label", main = "Line Plot Example")
