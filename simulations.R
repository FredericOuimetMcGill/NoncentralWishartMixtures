########################################################################################################
## On noncentral Wishart mixtures of noncentral Wisharts and their application to test random effects ##
## in factorial design models (Parallelized version)                                                  ##
########################################################################################################

## Written by Frederic Ouimet (February 2025)

path <- getwd()

############################
## Load required packages ##
############################

# These libraries will be loaded on both the master and the worker processes.
library("CholWishart")        # for the multivariate gamma function
library("cubature")           # for multidimensional integration (adaptIntegrate)
library("expm")               # for matrix square roots (if needed)
library("LaplacesDemon")      # for the Wishart and inverse Wishart distributions
library("MASS")               # for mvrnorm
library("matrixcalc")         # for is.positive.definite
library("matrixsampling")     # for the matrix-type-II Beta distribution
library("parallel")           # for parallel computing

######################################
## Determine number of cores to use ##
######################################

# Automatically detect the number of available cores and use them all.
num_cores <- detectCores()
cat("Using", num_cores, "cores for the simulation.\n")

# Create the cluster
cl <- makeCluster(num_cores)

####################################
## Export libraries and functions ##
####################################

# Load required libraries on each worker.
clusterEvalQ(cl, {
  library("CholWishart")
  library("cubature")
  library("expm")
  library("LaplacesDemon")
  library("MASS")
  library("matrixcalc")
  library("matrixsampling")
  NULL
})

##############################################
## Global simulation parameters and helpers ##
##############################################

# Set seed for reproducibility (only the master; workers get their own streams)
set.seed(12345)

# Model and simulation parameters:
d <- 2  # dimension of each observation (d >= 2)
a <- 3  # number of levels for factor A (a - 1 > d - 1)
b <- 3  # number of levels for factor B (b - 1 > d - 1)
n <- 5  # number of replicates per cell
nrep <- 10000  # Number of replications

# Overall mean (mu) and error covariance (Sigma)
mu <- rep(0, d)
Sigma <- diag(1, d)

######################
## Helper functions ##
######################

# Function to symmetrize a matrix X
symmetrize <- function(X) {
  return((X + t(X)) / 2)
}

# Function to compute the square root of a SPSD matrix P
matrix_sqrt <- function(P) {
  eigP <- eigen(P)
  sqrt_vals <- sqrt(eigP$values)
  return(eigP$vectors %*% diag(sqrt_vals) %*% t(eigP$vectors))
}

# Function to compute the conjugation R_P = P^(1/2) * R * P^(1/2)
conjugate <- function(R, P) {
  P_sqrt <- matrix_sqrt(P)
  return(P_sqrt %*% R %*% P_sqrt)
}

# Matrix-variate beta type II density function for a SPD matrix X
dmatrixbeta_typeII <- function(X, a, b) { # a, b > (d - 1)/2
  num <- mvgamma(a + b, d)
  den <- mvgamma(a, d) * mvgamma(b, d)
  det_X <- det(X)
  det_I_plus_X <- det(diag(d) + X)
  density <- num / den * det_X ^ (a - (d + 1) / 2) * det_I_plus_X ^ (-(a + b))
  return(density)
}

########################################################
## Test: Does the matrix beta density integrate to 1? ##
########################################################

## Coordinate functions for numerical integration

# Rotation matrix (2x2) given an angle theta.
rotation_matrix <- function(theta) {
  return(matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, byrow = TRUE))
}

# Construct a 2x2 symmetric matrix X from parameters theta, lambda1, lambda2.
construct_X <- function(theta, lambda1, lambda2) {
  R <- rotation_matrix(theta)
  diag_lambda <- diag(c(lambda1, lambda2))
  X <- R %*% diag_lambda %*% t(R)
  return(symmetrize(X)) # to ensure symmetry
}

# ## Test that the beta density integrates to 1
# 
# # Arbitrarily chosen parameters
# a_param <- 3
# b_param <- 4
# 
# # Define the integrand for the 3-dimensional integration
# integrand_test <- function(vars) {
#   theta <- vars[1]
#   lambda1 <- vars[2]
#   lambda2 <- vars[3]
# 
#   # Construct the symmetric positive-definite matrix X
#   X <- construct_X(theta, lambda1, lambda2)
# 
#   # Compute the density at X
#   fX <- dmatrixbeta_typeII(X, a_param, b_param)
# 
#   # Jacobian adjustment for the change of variables
#   jacobian <- abs(lambda1 - lambda2) / 4
# 
#   return(fX * jacobian)
# }
# 
# # Integration limits: theta in [0, 2*pi], lambda1 and lambda2 in (0, Inf)
# lowerLimits <- c(0, 0, 0)
# upperLimits <- c(2 * pi, Inf, Inf)
# 
# # Perform the integration
# integration_result <- adaptIntegrate(integrand_test, lowerLimit = lowerLimits, upperLimit = upperLimits, tol = 1e-3)
# 
# # Display the result
# cat("Integral of the matrix beta density (should be close to 1):", integration_result$integral, "\n")

########################################################################
## p-value computation for a 2x2 SPD test statistic (via Monte Carlo) ##
########################################################################

pvalue_matrix_beta_MC <- function(T_obs, a_param, b_param, n_MC = 10000) {
  f_obs <- dmatrixbeta_typeII(T_obs, a_param, b_param)
  
  SPD_list <- lapply(1:n_MC, function(x) {
    X <- matrixsampling::rmatrixbeta(1, d, a_param, b_param, def = 2, checkSymmetry = TRUE)[,,1]
    symmetrize(as.matrix(solve(diag(1, d) - X) %*% X))
  })
  
  indicators <- sapply(SPD_list, function(X) {
    fX <- dmatrixbeta_typeII(X, a_param, b_param)
    as.numeric(fX >= f_obs)
  })
  
  return(1 - mean(indicators))
}

###############################################
## Export variables and functions to workers ##
###############################################

clusterExport(cl, varlist = c("d", "a", "b", "n", "nrep", "mu", "Sigma",
                              "symmetrize", "matrix_sqrt", "conjugate",
                              "dmatrixbeta_typeII", "rotation_matrix",
                              "construct_X", "pvalue_matrix_beta_MC"))

#########################################################
## Simulation function: one replication (simulate_rep) ##
#########################################################

simulate_rep <- function(rep) {
  rep_start <- Sys.time()
  
  # Generate data:
  # Model: Y_{ijk} = mu + alpha_i + beta_j + (alpha*beta)_{ij} + epsilon_{ijk}
  Y_array <- array(NA, dim = c(a, b, n, d))
  for (i in 1:a) {
    alpha_i <- 0  # under H0^A
    for (j in 1:b) {
      beta_j <- 0  # under H0^B
      interaction_ij <- 0  # under H0^{AB}
      for (k in 1:n) {
        eps_ijk <- MASS::mvrnorm(1, mu = rep(0, d), Sigma = Sigma)
        Y_array[i, j, k, ] <- mu + alpha_i + beta_j + interaction_ij + eps_ijk
      }
    }
  }
  
  # Compute means and sums-of-cross-products (SCPs)
  all_obs <- matrix(Y_array, ncol = d)
  overall_mean <- colMeans(all_obs)
  
  row_means <- matrix(0, nrow = a, ncol = d)
  col_means <- matrix(0, nrow = b, ncol = d)
  
  for (i in 1:a) {
    row_means[i, ] <- colMeans(matrix(Y_array[i, , , ], ncol = d))
  }
  for (j in 1:b) {
    col_means[j, ] <- colMeans(matrix(Y_array[, j, , ], ncol = d))
  }
  
  # Initialize sum-of-cross-products matrices:
  SCP_A <- matrix(0, nrow = d, ncol = d)
  SCP_B <- matrix(0, nrow = d, ncol = d)
  SCP_AB <- matrix(0, nrow = d, ncol = d)
  SCP_E <- matrix(0, nrow = d, ncol = d)
  
  # SCP for Factor A:
  for (i in 1:a) {
    diff <- matrix(row_means[i, ] - overall_mean, ncol = 1)
    SCP_A <- SCP_A + (b * n) * (diff %*% t(diff))
  }
  
  # SCP for Factor B:
  for (j in 1:b) {
    diff <- matrix(col_means[j, ] - overall_mean, ncol = 1)
    SCP_B <- SCP_B + (a * n) * (diff %*% t(diff))
  }
  
  # SCP for Interaction AB:
  for (i in 1:a) {
    for (j in 1:b) {
      cell_ij <- matrix(Y_array[i, j, , ], ncol = d)
      cell_mean <- colMeans(cell_ij)
      diff <- matrix(cell_mean - row_means[i, ] - col_means[j, ] + overall_mean, ncol = 1)
      SCP_AB <- SCP_AB + n * (diff %*% t(diff))
    }
  }
  
  # SCP for Error:
  for (i in 1:a) {
    for (j in 1:b) {
      cell_ij <- matrix(Y_array[i, j, , ], ncol = d)
      cell_mean <- colMeans(cell_ij)
      for (k in 1:n) {
        diff <- matrix(cell_ij[k, ] - cell_mean, ncol = 1)
        SCP_E <- SCP_E + diff %*% t(diff)
      }
    }
  }
  
  # Compute test statistics:
  SCP_E_norm <- solve(conjugate(SCP_E, Sigma))
  
  T_A  <- conjugate(conjugate(SCP_A, Sigma), SCP_E_norm)
  T_B  <- conjugate(conjugate(SCP_B, Sigma), SCP_E_norm)
  T_AB <- conjugate(conjugate(SCP_AB, Sigma), SCP_E_norm)
  
  # Compute p-values via Monte Carlo:
  p_A  <- pvalue_matrix_beta_MC(T_A,  (a - 1) / 2, a * b * (n - 1) / 2)
  p_B  <- pvalue_matrix_beta_MC(T_B,  (b - 1) / 2, a * b * (n - 1) / 2)
  p_AB <- pvalue_matrix_beta_MC(T_AB, (a - 1) * (b - 1) / 2, a * b * (n - 1) / 2)
  
  rep_end <- Sys.time()
  rep_time <- as.numeric(difftime(rep_end, rep_start, units = "mins"))
  
  # Record p-values and rejection indicators for thresholds 0.01, 0.05, 0.10
  return(list(
    p_A = p_A,
    p_B = p_B,
    p_AB = p_AB,
    p_A_rate = c(as.numeric(p_A < 0.01), as.numeric(p_A < 0.05), as.numeric(p_A < 0.10)),
    p_B_rate = c(as.numeric(p_B < 0.01), as.numeric(p_B < 0.05), as.numeric(p_B < 0.10)),
    p_AB_rate = c(as.numeric(p_AB < 0.01), as.numeric(p_AB < 0.05), as.numeric(p_AB < 0.10)),
    rep_time = rep_time
  ))
}

########################################################
## Run the simulation in parallel (over replications) ##
########################################################

total_start_time <- Sys.time()

# Run the replications in parallel
results_list <- parLapply(cl, 1:nrep, simulate_rep)

total_end_time <- Sys.time()
total_time <- as.numeric(difftime(total_end_time, total_start_time, units = "mins"))

# Stop the cluster
stopCluster(cl)

######################################
## Combine and save simulation data ##
######################################

# Extract rates from the results list:
results_A  <- do.call(rbind, lapply(results_list, function(x) x$p_A_rate))
results_B  <- do.call(rbind, lapply(results_list, function(x) x$p_B_rate))
results_AB <- do.call(rbind, lapply(results_list, function(x) x$p_AB_rate))
rep_times  <- sapply(results_list, function(x) x$rep_time)

# Create raw data frames for each test:
raw_A <- data.frame(Test = rep("Factor A", nrep), results_A)
raw_B <- data.frame(Test = rep("Factor B", nrep), results_B)
raw_AB <- data.frame(Test = rep("Interaction AB", nrep), results_AB)
raw_table <- rbind(raw_A, raw_B, raw_AB)

# Save raw rejection indicators to CSV:
write.csv(raw_table, file = "raw_table.csv", row.names = FALSE)

########################################
## Summarize and save overall results ##
########################################

# Function to summarize (mean and sd) the 0/1 indicators.
summarize_results <- function(mat, testname) {
  res <- data.frame(
    Test = testname,
    Threshold = c("0.01", "0.05", "0.10"),
    Mean = apply(mat, 2, mean, na.rm = TRUE),
    SE = apply(mat, 2, sd, na.rm = TRUE) / sqrt(nrow(mat))
  )
  res
}

res_A  <- summarize_results(results_A, "Factor A")
res_B  <- summarize_results(results_B, "Factor B")
res_AB <- summarize_results(results_AB, "Interaction AB")
final_table <- rbind(res_A, res_B, res_AB)

cat("\nSummary table of rejection rates (proportion of p-values below threshold):\n")
print(final_table)

# Save summary table to CSV:
write.csv(final_table, file = "summary_table.csv", row.names = FALSE)

cat("\nAverage replication time (mins):", mean(rep_times), "\n")
cat("Total simulation time (mins):", total_time, "\n")

###############################
## Uniform tests on p_values ##
###############################

# Load the nortest package for Anderson-Darling test
library("nortest")  # For the Anderson-Darling test

# Extract p-values from the results list:
p_A_values  <- do.call(rbind, lapply(results_list, function(x) x$p_A))
p_B_values  <- do.call(rbind, lapply(results_list, function(x) x$p_B))
p_AB_values <- do.call(rbind, lapply(results_list, function(x) x$p_AB))

##################################
## Save raw p-values to CSV file ##
##################################

# Ensure p-values are stored as vectors
p_A_values <- as.vector(p_A_values)
p_B_values <- as.vector(p_B_values)
p_AB_values <- as.vector(p_AB_values)

# Create raw data frames for each test:
raw_p_A  <- data.frame(Test = rep("Factor A", length(p_A_values)), p_value = p_A_values)
raw_p_B  <- data.frame(Test = rep("Factor B", length(p_B_values)), p_value = p_B_values)
raw_p_AB <- data.frame(Test = rep("Interaction AB", length(p_AB_values)), p_value = p_AB_values)

# Combine into one table
raw_p_table <- rbind(raw_p_A, raw_p_B, raw_p_AB)

# Save raw p-values to CSV:
write.csv(raw_p_table, file = "raw_p_values.csv", row.names = FALSE)

############################################
## Perform uniformity tests on p-values  ##
############################################

# Kolmogorov-Smirnov (K-S) Test for Uniform(0,1)
ks_test_A <- ks.test(p_A_values, "punif", 0, 1)
ks_test_B <- ks.test(p_B_values, "punif", 0, 1)
ks_test_AB <- ks.test(p_AB_values, "punif", 0, 1)

############################################
## Save summary of uniformity test results ##
############################################

# Function to format test results
summarize_uniformity_test <- function(ks_test, testname) {
  data.frame(
    Test = testname,
    KS_Statistic = ks_test$statistic,
    KS_p_value = ks_test$p.value
  )
}

# Create summary table
summary_uniformity_table <- rbind(
  summarize_uniformity_test(ks_test_A, "Factor A"),
  summarize_uniformity_test(ks_test_B, "Factor B"),
  summarize_uniformity_test(ks_test_AB, "Interaction AB")
)

# Print summary
cat("\nSummary table of uniformity test results:\n")
print(summary_uniformity_table)

# Save summary table to CSV:
write.csv(summary_uniformity_table, file = "summary_uniformity_tests.csv", row.names = FALSE)

######################################
## Print test results to console ##
######################################

cat("\nK-S Test Results:\n")
print(ks_test_A)
print(ks_test_B)
print(ks_test_AB)


############################################################
############################################################
##--------------------------------------------------------##
##    Univariate (d=1) Two-way ANOVA Simulation Code      ##
##--------------------------------------------------------##
############################################################
############################################################

set.seed(12345)  # for reproducibility

# Dimensions and design:
d <- 1   # dimension of each observation (now a scalar)
a <- 3   # number of levels for factor A
b <- 3   # number of levels for factor B
n <- 5   # number of replicates per cell

# Number of Monte Carlo replications
nrep <- 100000

# True grand mean:
mu <- 0

# Error variance:
sigma_sq <- 1

# Storage for rejection indicators:
# For each test (A, B, AB) and for each threshold, we record 1 if p-value < threshold, else 0.
results <- list(
  A  = matrix(0, nrow = nrep, ncol = 3),
  B  = matrix(0, nrow = nrep, ncol = 3),
  AB = matrix(0, nrow = nrep, ncol = 3)
)
colnames(results$A) <- colnames(results$B) <- colnames(results$AB) <- c("p<0.01", "p<0.05", "p<0.10")

# To record the time (in minutes) of each replication:
rep_times <- numeric(nrep)

#################################
## Main Monte Carlo simulation ##
#################################

total_start_time <- Sys.time()

for (rep in 1:nrep) {
  rep_start <- Sys.time()
  
  # ------------------------------------------------------------------------------
  # 1) Generate data:
  #    Model: Y_{ijk} = mu + alpha_i + beta_j + (alpha*beta)_{ij} + epsilon_{ijk},
  #    where d=1 => all of these are scalars, and under H0 all factor effects = 0.
  # ------------------------------------------------------------------------------
  
  # Array to store Y_{ijk}: dimension (a, b, n).
  Y_array <- array(NA, dim = c(a, b, n))
  
  for (i in 1:a) {
    for (j in 1:b) {
      for (k in 1:n) {
        # Under H0, alpha_i, beta_j, and interaction_ij are all zero.
        # So Y_{ijk} = mu + random error
        Y_array[i, j, k] <- mu + rnorm(1, mean = 0, sd = sqrt(sigma_sq))
      }
    }
  }
  
  # ------------------------------------------------------------------------
  # 2) Compute various means and sums of squares for the ANOVA decomposition
  # ------------------------------------------------------------------------
  
  # Overall mean
  all_obs <- as.vector(Y_array)  # unroll everything into a vector
  overall_mean <- mean(all_obs)
  
  # Row means (factor A)
  row_means <- numeric(a)
  for (i in 1:a) {
    row_means[i] <- mean(Y_array[i, , ])
  }
  
  # Column means (factor B)
  col_means <- numeric(b)
  for (j in 1:b) {
    col_means[j] <- mean(Y_array[, j, ])
  }
  
  # Cell means for interaction
  cell_means <- matrix(0, nrow = a, ncol = b)
  for (i in 1:a) {
    for (j in 1:b) {
      cell_means[i, j] <- mean(Y_array[i, j, ])
    }
  }
  
  # ------------------------------------
  # 3) Compute means and sums-of-squares
  # ------------------------------------
  
  # Compute SS_A
  SS_A <- b * n * sum((row_means - overall_mean)^2)
  
  # Compute SS_B
  SS_B <- a * n * sum((col_means - overall_mean)^2)
  
  # Compute SS_AB
  SS_AB <- 0
  for (i in 1:a) {
    for (j in 1:b) {
      SS_AB <- SS_AB + n * (cell_means[i, j] - row_means[i] - col_means[j] + overall_mean)^2
    }
  }
  
  # Compute SS_E
  SS_E <- 0
  for (i in 1:a) {
    for (j in 1:b) {
      for (k in 1:n) {
        SS_E <- SS_E + (Y_array[i, j, k] - cell_means[i, j])^2
      }
    }
  }
  
  # ----------------------------------------------
  # 4) Compute classical F-statistics and p-values
  # ----------------------------------------------
  
  # Degrees of freedom
  df_A  <- a - 1
  df_B  <- b - 1
  df_AB <- (a - 1) * (b - 1)
  df_E  <- a * b * (n - 1)
  
  # Mean squares
  MS_A  <- SS_A  / df_A
  MS_B  <- SS_B  / df_B
  MS_AB <- SS_AB / df_AB
  MS_E  <- SS_E  / df_E
  
  # F-statistics
  F_A  <- MS_A  / MS_E
  F_B  <- MS_B  / MS_E
  F_AB <- MS_AB / MS_E
  
  # p-values from the F distribution
  p_A  <- 1 - pf(F_A,  df_A,  df_E)
  p_B  <- 1 - pf(F_B,  df_B,  df_E)
  p_AB <- 1 - pf(F_AB, df_AB, df_E)
  
  # -----------------------------------------------------------------
  # 5) Record whether p-values fall below 0.01, 0.05, 0.10 thresholds
  # -----------------------------------------------------------------
  results$A[rep, ]  <- c(as.numeric(p_A < 0.01),  as.numeric(p_A < 0.05),  as.numeric(p_A < 0.10))
  results$B[rep, ]  <- c(as.numeric(p_B < 0.01),  as.numeric(p_B < 0.05),  as.numeric(p_B < 0.10))
  results$AB[rep, ] <- c(as.numeric(p_AB < 0.01), as.numeric(p_AB < 0.05), as.numeric(p_AB < 0.10))
  
  # 6) Record replication time
  rep_end <- Sys.time()
  rep_times[rep] <- as.numeric(difftime(rep_end, rep_start, units = "mins"))
  
  if (rep %% 50 == 0) {
    cat("Completed replication", rep, "of", nrep, "\n")
  }
}

total_end_time <- Sys.time()
total_time <- as.numeric(difftime(total_end_time, total_start_time, units = "mins"))
mean_time  <- mean(rep_times)

###################################
## Save raw rejection indicators ##
###################################

# Create a raw data frame for each test (adding a Test column) without summarizing
raw_A  <- data.frame(Test = rep("Factor A",        nrep), results$A)
raw_B  <- data.frame(Test = rep("Factor B",        nrep), results$B)
raw_AB <- data.frame(Test = rep("Interaction AB",  nrep), results$AB)

# Combine into one raw table
raw_table <- rbind(raw_A, raw_B, raw_AB)

# Save the raw table to a CSV file:
write.csv(raw_table, file = "raw_table_d1.csv", row.names = FALSE)

################################
## Summarize and save results ##
################################

# Function to summarize (mean and standard error) the 0/1 indicators
summarize_results <- function(mat, testname) {
  nrep <- nrow(mat)
  data.frame(
    Test      = testname,
    Threshold = c("0.01", "0.05", "0.10"),
    Mean      = apply(mat, 2, mean),
    SE        = apply(mat, 2, sd) / sqrt(nrep)
  )
}

res_A  <- summarize_results(results$A,  "Factor A")
res_B  <- summarize_results(results$B,  "Factor B")
res_AB <- summarize_results(results$AB, "Interaction AB")
final_table <- rbind(res_A, res_B, res_AB)

cat("\nSummary table of rejection rates:\n")
print(final_table)

# Save the summary table to a CSV file:
write.csv(final_table, file = "summary_table_d1.csv", row.names = FALSE)

##############################
## Print timing information ##
##############################

cat("\nTiming Information:\n")
cat("Mean time per replication (minutes):", round(mean_time, 4), "\n")
cat("Total time (minutes):", round(total_time, 4), "\n")
