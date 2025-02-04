###################################################################################################################################
## On noncentral Wishart mixtures of noncentral Wisharts and their application to test random effects in factorial design models ##
###################################################################################################################################

## Written by Frederic Ouimet (February 2025)

##################
## Set the path ##
##################

path <- file.path(
  "C://Users//fred1//Dropbox//Ouimet_Genest_projects//GMOR_2025_noncentral_Wishart_mixture//_simulations",
  fsep = .Platform$file.sep
)
# path <- getwd()
setwd(path)

###########################################
## Set-up: load libraries and parameters ##
###########################################

# Load required packages

library("CholWishart")        # for the multivariate gamma function
library("cubature")           # for multidimensional integration (adaptIntegrate)
library("expm")               # for matrix square roots (if needed)
library("LaplacesDemon")      # for the Wishart and inverse Wishart distributions
library("MASS")               # for mvrnorm
library("matrixcalc")         # for is.positive.definite

# -------------------------------
# Model and simulation parameters
# -------------------------------

set.seed(12345)  # for reproducibility

# Dimensions and design:
d <- 2               # dimension of each observation (2-dimensional)
a <- 3               # number of levels for factor A (a - 1 > d - 1)
b <- 3               # number of levels for factor B (b - 1 > d - 1)
n <- 5               # number of replicates per cell

# Overall mean (2-dimensional):
mu <- c(1, 2)

# Error covariance matrix (Sigma) (2x2):
Sigma <- matrix(c(1, 0.3,
                  0.3, 1), nrow = 2)

# Covariance matrices of the factors (2x2):
Sigma_alpha     <- matrix(c(1, -0.5, -0.5, 1), nrow = 2, ncol = 2)  # for factor A
Sigma_beta      <- matrix(c(1, 0.2, 0.2, 1), nrow = 2, ncol = 2)  # for factor B
Sigma_alpha_beta<- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)  # for the interaction

# Integration parameters for p-value calculation:
# We will parameterize a 2x2 SPD matrix X by:
# X = R(theta) %*% diag(lambda1, lambda2) %*% R(theta)^T,
# with theta in [0, 2*pi] and lambda1, lambda2 in [0, Lmax].
Lmax <- 1      # upper limit for eigenvalues in integration
tol_int <- 1e-3 # tolerance for the integration (adjust for speed/accuracy)

# Replication parameters:
nrep <- 1000  # number of replications

# Storage for rejection indicators:
# For each test (A, B, AB) and for each threshold, we record 1 if p-value < threshold, else 0.
results <- list(
  A = matrix(0, nrow = nrep, ncol = 3),
  B = matrix(0, nrow = nrep, ncol = 3),
  AB = matrix(0, nrow = nrep, ncol = 3)
)
colnames(results$A) <- colnames(results$B) <- colnames(results$AB) <- c("p<0.01", "p<0.05", "p<0.10")

# To record the time (in minutes) of each replication:
rep_times <- numeric(nrep)

################################################################################
## Helper functions: conjugation, matrix square-root inverse and beta density ##
################################################################################

# Function to compute the square root of a symmetric positive–semidefinite matrix P.
matrix_sqrt <- function(P) {
  # Check if P is symmetric
  if (!isSymmetric(P)) {
    stop("Matrix P must be symmetric.")
  }
  # Compute eigen-decomposition: P = Q Λ Qᵀ
  eigenP <- eigen(P)
  # Ensure nonnegative eigenvalues (in case of numerical issues)
  sqrt_eigen <- sqrt(pmax(eigenP$values, 0))
  # Compute the square root: P^(1/2) = Q * sqrt(Λ) * Qᵀ
  P_sqrt <- eigenP$vectors %*% diag(sqrt_eigen) %*% t(eigenP$vectors)
  return(P_sqrt)
}

# Function to compute the conjugation R_P = P^(1/2) * R * P^(1/2)
conjugate_matrix <- function(R, P) {
  # Compute the square root of P
  P_sqrt <- matrix_sqrt(P)
  # Return the conjugated matrix
  R_P <- P_sqrt %*% R %*% P_sqrt
  return(R_P)
}

# Function to compute the inverse square root of a SPD matrix Sigma:
inv_sqrt <- function(Sigma) {
  eig <- eigen(Sigma)
  # Avoid division by zero:
  inv_sqrt_vals <- 1/sqrt(eig$values)
  return(eig$vectors %*% diag(inv_sqrt_vals) %*% t(eig$vectors))
}

# Matrix-variate beta type II density function for a 2x2 SPD matrix X.
# X is 2x2, and parameters a and b must be > (d-1)/2 = 0.5.
dmatrixbeta_typeII <- function(X, a, b) {
  d <- nrow(X)
  num <- CholWishart::mvgamma(a + b, d)
  den <- CholWishart::mvgamma(a, d) * CholWishart::mvgamma(b, d)
  detX <- det(X)
  detI_plus_X <- det(diag(d) + X)
  density <- (num/den) * (detX)^(a - (d + 1)/2) * (detI_plus_X)^(-(a + b))
  return(density)
}

##########################################
## Coordinate functions for integration ##
##########################################

# Rotation matrix (2x2) given an angle theta.
rotation_matrix <- function(theta) {
  return(matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, byrow = TRUE))
}

# Construct a 2x2 symmetric matrix X from parameters theta, lambda1, lambda2.
construct_X <- function(theta, lambda1, lambda2) {
  R <- rotation_matrix(theta)
  diag_lambda <- diag(c(lambda1, lambda2))
  X <- R %*% diag_lambda %*% t(R)
  X <- (X + t(X)) / 2 # to ensure symmetry
  return(X)
}

######################################################
## p-value computation for a 2x2 SPD test statistic ##
######################################################

# Given an observed 2x2 SPD test statistic T_obs, and beta parameters (a_param, b_param),
# compute the p-value as 1 minus the integrated probability of all matrices X for which
# dmatrixbeta_typeII(X, a_param, b_param) >= dmatrixbeta_typeII(T_obs, a_param, b_param).
pvalue_matrix_beta <- function(T_obs, a_param, b_param) {
  # Evaluate the density at the observed test statistic:
  f_obs <- dmatrixbeta_typeII(T_obs, a_param, b_param)
  
  # Define the integrand for the 3-dimensional integration.
  # We parameterize an SPD matrix X by (theta, lambda1, lambda2).
  integrand <- function(vars) {
    theta <- vars[1]
    lambda1 <- vars[2]
    lambda2 <- vars[3]
    # Construct X:
    X <- construct_X(theta, lambda1, lambda2)
    # Compute the density at X:
    fX <- dmatrixbeta_typeII(X, a_param, b_param)
    # Indicator: 1 if fX >= f_obs, 0 otherwise.
    indicator <- as.numeric(fX >= f_obs)
    # Jacobian adjustment for the change of variables:
    jacobian <- abs(lambda1 - lambda2) / 4
    return(fX * indicator * jacobian)
  }
  
  # Integration limits: theta in [0,2*pi], lambda1 and lambda2 in [0, Lmax]
  lowerLimits <- c(0, 0, 0)
  upperLimits <- c(2*pi, Lmax, Lmax)
  
  # Perform the integration using adaptIntegrate (cubature package)
  res <- try(adaptIntegrate(integrand, lowerLimit = lowerLimits, upperLimit = upperLimits, tol = tol_int),
             silent = TRUE)
  if(inherits(res, "try-error")) {
    integrated_value <- NA
  } else {
    integrated_value <- res$integral
  }
  # p-value is one minus the integrated probability:
  p_val <- 1 - integrated_value
  # In case numerical integration returns p_val slightly < 0 or > 1, adjust:
  p_val <- max(min(p_val, 1), 0)
  return(p_val)
}

##########################
## Main simulation loop ##
##########################

total_start_time <- Sys.time()

for (rep in 1:nrep) {
  rep_start <- Sys.time()
  
  # --------------------------------------------------------------------------
  # Generate data:
  # Model: Y_{ijk} = mu + alpha_i + beta_j + (alpha*beta)_{ij} + epsilon_{ijk}
  # For testing under H0, set alpha_i, beta_j, and interaction to 0.
  # Each Y_{ijk} is 2-dimensional.
  # --------------------------------------------------------------------------
  
  # Preallocate a 4D array: dimensions (a, b, n, d)
  Y_array <- array(NA, dim = c(a, b, n, d))
  
  # For each level of factor A and B, generate n observations.
  for (i in 1:a) {
    # Under H0, alpha_i = 0 (2-dim zero vector)
    alpha_i <- rep(0, d)
    for (j in 1:b) {
      # Under H0, beta_j = 0 and interaction = 0.
      beta_j <- rep(0, d)
      interaction_ij <- rep(0, d)
      # For each replicate in cell (i,j), generate a 2-dim observation:
      for (k in 1:n) {
        # epsilon ~ N(0, Sigma)
        eps <- mvrnorm(1, mu = rep(0, d), Sigma = Sigma)
        Y_array[i, j, k, ] <- mu + alpha_i + beta_j + interaction_ij + eps
      }
    }
  }
  
  # -----------------------------------------
  # Compute means and sums-of-cross-products:
  # -----------------------------------------
  
  # Compute overall mean (2-dimensional vector)
  all_obs <- matrix(Y_array, ncol = d)  # each row is an observation
  overall_mean <- colMeans(all_obs)
  
  # Initialize SCP matrices as 2x2 zero matrices:
  SCP_A <- matrix(0, nrow = d, ncol = d)
  SCP_B <- matrix(0, nrow = d, ncol = d)
  SCP_AB <- matrix(0, nrow = d, ncol = d)
  SCP_E <- matrix(0, nrow = d, ncol = d)
  
  # Compute cell means (for each (i,j)) and store row and column means:
  cell_means <- matrix(NA, nrow = a * b, ncol = d)  # each row is a cell mean
  cell_index <- matrix(NA, nrow = a, ncol = b)      # to index into cell_means
  counter <- 1
  for (i in 1:a) {
    for (j in 1:b) {
      cell_ij <- matrix(Y_array[i, j, , ], ncol = d)
      cell_mean_ij <- colMeans(cell_ij)
      cell_means[counter, ] <- cell_mean_ij
      cell_index[i, j] <- counter
      counter <- counter + 1
    }
  }
  
  # Row means for factor A (average over cells in row i):
  row_means <- matrix(NA, nrow = a, ncol = d)
  for (i in 1:a) {
    # Average over all j and replicates:
    temp <- NULL
    for (j in 1:b) {
      cell_ij <- matrix(Y_array[i, j, , ], ncol = d)
      temp <- rbind(temp, cell_ij)
    }
    row_means[i, ] <- colMeans(temp)
  }
  
  # Column means for factor B:
  col_means <- matrix(NA, nrow = b, ncol = d)
  for (j in 1:b) {
    temp <- NULL
    for (i in 1:a) {
      cell_ij <- matrix(Y_array[i, j, , ], ncol = d)
      temp <- rbind(temp, cell_ij)
    }
    col_means[j, ] <- colMeans(temp)
  }
  
  # Compute SCP for factor A:
  for (i in 1:a) {
    diff <- matrix(row_means[i, ] - overall_mean, ncol = 1)
    SCP_A <- SCP_A + (b * n) * (diff %*% t(diff))
  }
  
  # Compute SCP for factor B:
  for (j in 1:b) {
    diff <- matrix(col_means[j, ] - overall_mean, ncol = 1)
    SCP_B <- SCP_B + (a * n) * (diff %*% t(diff))
  }
  
  # Compute SCP for the interaction:
  for (i in 1:a) {
    for (j in 1:b) {
      cell_ij <- matrix(Y_array[i, j, , ], ncol = d)
      cell_mean <- colMeans(cell_ij)
      diff <- matrix(cell_mean - row_means[i, ] - col_means[j, ] + overall_mean, ncol = 1)
      SCP_AB <- SCP_AB + n * (diff %*% t(diff))
    }
  }
  
  # Compute error SCP:
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
  
  # ------------------------
  # Compute test statistics:
  # ------------------------
  
  SCP_E_norm <- inv_sqrt(conjugate_matrix(SCP_E, Sigma))
  
  T_A  <- SCP_E_norm * conjugate_matrix(SCP_A, Sigma) * SCP_E_norm
  T_B  <- SCP_E_norm * conjugate_matrix(SCP_B, Sigma) * SCP_E_norm
  T_AB <- SCP_E_norm * conjugate_matrix(SCP_AB, Sigma) * SCP_E_norm
  
  # ---------------------------------------------------------
  # Compute p-values for each test via numerical integration.
  # ---------------------------------------------------------
  
  p_A  <- pvalue_matrix_beta(T_A, (a - 1) / 2,  a * b * (n - 1) / 2)
  p_B  <- pvalue_matrix_beta(T_B, (b - 1) / 2,  a * b * (n - 1) / 2)
  p_AB <- pvalue_matrix_beta(T_AB, (a - 1) * (b - 1) / 2, a * b * (n - 1) / 2)
  
  # -------------------------------------------------------------------
  # Record whether p-values fall below the thresholds 0.01, 0.05, 0.10.
  # -------------------------------------------------------------------
  results$A[rep, ]  <- c(as.numeric(p_A < 0.01), as.numeric(p_A < 0.05), as.numeric(p_A < 0.10))
  results$B[rep, ]  <- c(as.numeric(p_B < 0.01), as.numeric(p_B < 0.05), as.numeric(p_B < 0.10))
  results$AB[rep, ] <- c(as.numeric(p_AB < 0.01), as.numeric(p_AB < 0.05), as.numeric(p_AB < 0.10))
  
  # -------------------------------------
  # Record replication time (in minutes).
  # -------------------------------------
  rep_end <- Sys.time()
  rep_times[rep] <- as.numeric(difftime(rep_end, rep_start, units = "mins"))
  
  if(rep %% 50 == 0) {
    cat("Completed replication", rep, "of", nrep, "\n")
  }
}

total_end_time <- Sys.time()
total_time <- as.numeric(difftime(total_end_time, total_start_time, units = "mins"))
mean_time  <- mean(rep_times)

################################
## Summarize and save results ##
################################

# Function to summarize (mean and sd) the 0/1 indicators.
summarize_results <- function(mat, testname) {
  res <- data.frame(Test = testname,
                    Threshold = c("0.01", "0.05", "0.10"),
                    Mean = apply(mat, 2, mean, na.rm = TRUE),
                    SD = apply(mat, 2, sd, na.rm = TRUE))
  return(res)
}

res_A  <- summarize_results(results$A, "Factor A")
res_B  <- summarize_results(results$B, "Factor B")
res_AB <- summarize_results(results$AB, "Interaction AB")
final_table <- rbind(res_A, res_B, res_AB)

cat("\nSummary table of rejection rates (proportion of p-values below threshold):\n")
print(final_table)

# Save the summary table to a CSV file:
write.csv(final_table, file = "rejection_rates_summary_multivariate.csv", row.names = FALSE)

##############################
## Print timing information ##
##############################

cat("\nTiming Information:\n")
cat("Mean time per replication (minutes):", round(mean_time, 4), "\n")
cat("Total time (minutes):", round(total_time, 4), "\n")
