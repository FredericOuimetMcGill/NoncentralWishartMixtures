#####################################################################
## Exact Beta-II MANOVA vs three separate two-way ANOVAs on NHANES
## Factors: Education (5 levels)  ×  MaritalStatus (6 levels)
## Responses: BMI, TotChol (d = 2)
#####################################################################

## Written by Frédéric Ouimet (July 2025)

## ----- 0. Parameters -------------------------------------------------------

set.seed(2025)                             # reproducibility
d <- 2                                     # dimension of the response
a <- 5 ; b <- 6 ; n <- 5                   # 5 × 6 × 5 balanced design
deci <- 8                                  # number of decimals for the p-values

## ----- 1. Load libraries ---------------------------------------------------

library(NHANES)                            # NHANES 2009–2012 micro-data
library(dplyr)                             # data-manipulation verbs
library(CholWishart)                       # mvgamma() for Beta-II density
library(matrixsampling)                    # rmatrixbeta() for Monte-Carlo draws

## ----- 2. Helper functions -------------------------------------------------

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

# Monte-Carlo p-value for a Beta II test statistic
pvalue_matrix_beta_MC <- function(T_obs, a_param, b_param, n_MC = 1e4) {
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

## ----- 3. Build balanced 5×6×n dataset ------------------------------------

nh_bal <- NHANES |>
  # keep rows with complete data in factors + responses
  filter(!is.na(Education), !is.na(MaritalStatus),
         !is.na(BMI), !is.na(TotChol)) |>
  # recode factors with the *exact* short labels used by the package
  mutate(
    Education = factor(trimws(Education), levels = c(
      "8th Grade",
      "9 - 11th Grade",
      "High School",
      "Some College",
      "College Grad"
    )),
    MaritalStatus = factor(trimws(MaritalStatus), levels = c(
      "Married",
      "LivePartner",
      "Divorced",
      "Separated",
      "Widowed",
      "NeverMarried"
    ))
  ) |>
  # form the 5 × 6 cells
  group_by(Education, MaritalStatus) |>
  filter(n() >= n) |>
  slice_head(n = n) |>
  ungroup()

stopifnot(
  n_distinct(nh_bal$Education)     == a,
  n_distinct(nh_bal$MaritalStatus) == b,
  nrow(nh_bal)                     == a * b * n   # 5 × 6 × 5 = 150 rows
)

# build response matrix Y (150 × 2)
Y <- as.matrix(nh_bal[, c("BMI", "TotChol")])
rownames(Y) <- with(
  nh_bal,
  paste(Education, MaritalStatus,
        ave(seq_len(nrow(nh_bal)), Education, MaritalStatus, FUN = seq_along),
        sep = "_")
)

## ----- 4. Compute MANOVA sums-of-outer-products ---------------------------

grand_mean <- colMeans(Y)                   # overall mean (length 2)
resp       <- c("BMI", "TotChol")           # response names

# means by Education (factor A)
row_means <- nh_bal |>
  group_by(Education) |>
  summarise(across(all_of(resp), mean), .groups = "drop") |>
  select(all_of(resp)) |> as.matrix()

# means by MaritalStatus (factor B)
col_means <- nh_bal |>
  group_by(MaritalStatus) |>
  summarise(across(all_of(resp), mean), .groups = "drop") |>
  select(all_of(resp)) |> as.matrix()

# means for each Education × MaritalStatus cell
cell_means <- nh_bal |>
  group_by(Education, MaritalStatus) |>
  summarise(across(all_of(resp), mean), .groups = "drop") |>
  arrange(Education, MaritalStatus) |>
  select(all_of(resp)) |> as.matrix()
dim(cell_means) <- c(a, b, d)               # 5 × 6 × 2 array

# initialise SOP matrices
SOP_A <- SOP_B <- SOP_AB <- SOP_E <- matrix(0, d, d)

# between-rows SOP (Education)
for(i in seq_len(a)){
  diff  <- matrix(row_means[i, ] - grand_mean, ncol = 1)
  SOP_A <- SOP_A + (b * n) * diff %*% t(diff)
}

# between-columns SOP (MaritalStatus)
for(j in seq_len(b)){
  diff  <- matrix(col_means[j, ] - grand_mean, ncol = 1)
  SOP_B <- SOP_B + (a * n) * diff %*% t(diff)
}

# interaction SOP (Education × MaritalStatus)
for(i in seq_len(a)) for(j in seq_len(b)){
  diff   <- matrix(cell_means[i, j, ] - row_means[i, ] -
                     col_means[j, ]     + grand_mean, ncol = 1)
  SOP_AB <- SOP_AB + n * diff %*% t(diff)
}

# residual / error SOP
for(r in seq_len(nrow(nh_bal))){
  i <- match(nh_bal$Education[r], levels(nh_bal$Education))
  j <- match(nh_bal$MaritalStatus[r], levels(nh_bal$MaritalStatus))
  diff <- matrix(Y[r, ] - cell_means[i, j, ], ncol = 1)
  SOP_E <- SOP_E + diff %*% t(diff)
}

## ----- 5. Form Beta-II test statistics ------------------------------------

Sigma_hat <- diag(d)
Vinv <- solve(conjugate(SOP_E, Sigma_hat))

# Test statistics
T_A   <- conjugate(conjugate(SOP_A, Sigma), Vinv)
T_B   <- conjugate(conjugate(SOP_B, Sigma), Vinv)
T_AB  <- conjugate(conjugate(SOP_AB, Sigma), Vinv)

# Monte-Carlo p-values
p_A  <- pvalue_matrix_beta_MC(T_A, (a - 1)/2, a * b * (n - 1)/2)
p_B  <- pvalue_matrix_beta_MC(T_B, (b - 1)/2, a * b * (n - 1)/2)
p_AB <- pvalue_matrix_beta_MC(T_AB, (a - 1)*(b - 1)/2, a * b * (n - 1)/2)

## ----- 6. Univariate two-way ANOVAs ---------------------------------------

get_p <- function(var, term){
  smry <- summary(aov(reformulate(c("Education", "MaritalStatus",
                                    "Education:MaritalStatus"), var),
                      data = nh_bal))
  smry[[1]][["Pr(>F)"]][term]   # 1 = Education, 2 = MaritalStatus, 3 = interaction
}

p_uni <- data.frame(
  Variable = resp,
  p_A      = sapply(resp, get_p, term = 1),
  p_B      = sapply(resp, get_p, term = 2),
  p_AB     = sapply(resp, get_p, term = 3),
  row.names = NULL
) |>
  mutate(across(where(is.numeric), \(x) round(x, deci)))

# ## ----- 6. Univariate two-way ANOVAs (from scratch, this is a test) --------
# 
# # helper: balanced-design two-way ANOVA p-values without aov()
# two_way_anova_p <- function(y, facA, facB, a, b, n) {
# 
#   y      <- as.numeric(y)
#   facA   <- factor(facA)        # Education  (a = 5 levels)
#   facB   <- factor(facB)        # Marital    (b = 6 levels)
# 
#   grand_mean <- mean(y)
# 
#   # cell means, row (A) means, column (B) means
#   mean_AB <- tapply(y, list(facA, facB), mean)
#   mean_A  <- rowMeans(mean_AB)
#   mean_B  <- colMeans(mean_AB)
# 
#   # sum-of-squares
#   SS_A  <- b * n * sum( (mean_A  - grand_mean)^2 )
#   SS_B  <- a * n * sum( (mean_B  - grand_mean)^2 )
#   SS_AB <- n     * sum( (mean_AB - outer(mean_A, mean_B, "+") + grand_mean)^2 )
#   SS_E  <- sum( ( y - mean_AB[cbind(as.numeric(facA), as.numeric(facB))] )^2 )
# 
#   # degrees of freedom
#   df_A  <- a - 1
#   df_B  <- b - 1
#   df_AB <- (a - 1) * (b - 1)
#   df_E  <- a * b * (n - 1)
# 
#   # mean squares and F statistics
#   MS_A  <- SS_A  / df_A
#   MS_B  <- SS_B  / df_B
#   MS_AB <- SS_AB / df_AB
#   MS_E  <- SS_E  / df_E
# 
#   F_A   <- MS_A  / MS_E
#   F_B   <- MS_B  / MS_E
#   F_AB  <- MS_AB / MS_E
# 
#   # p-values
#   p_A   <- pf(F_A, df_A, df_E, lower.tail = FALSE)
#   p_B   <- pf(F_B, df_B, df_E, lower.tail = FALSE)
#   p_AB  <- pf(F_AB, df_AB, df_E, lower.tail = FALSE)
# 
#   c(p_A, p_B, p_AB)
# }
# 
# # compute p-values for each response
# p_vals <- lapply(resp, \(v)
#                  two_way_anova_p(nh_bal[[v]],
#                                  nh_bal$Education,
#                                  nh_bal$MaritalStatus,
#                                  a, b, n) )
# 
# p_uni_from_scratch <- data.frame(
#   Variable = resp,
#   p_A  = sapply(p_vals, `[[`, 1),
#   p_B  = sapply(p_vals, `[[`, 2),
#   p_AB = sapply(p_vals, `[[`, 3),
#   row.names = NULL
# ) |>
#   mutate(across(where(is.numeric), \(x) round(x, deci)))

## ----- 7. Display results --------------------------------------------------

cat("\nExact Beta-II MANOVA p-values:\n")
print(formatC(c(p_A = p_A, p_B = p_B, p_AB = p_AB),
              format = "f", digits = deci))

cat("\nUnivariate two-way ANOVA p-values:\n")
print(p_uni)

# cat("\nUnivariate two-way ANOVA p-values (from scratch):\n")
# print(p_uni_from_scratch)
