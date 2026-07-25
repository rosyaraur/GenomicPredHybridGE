# ==============================================================================
# C++ ACCELERATED HYBRID PREDICTION PIPELINE
# Eigen-Decomposition RKHS + Custom RcppArmadillo Gibbs Sampler
# ==============================================================================

if (!require("Rcpp")) install.packages("Rcpp")
if (!require("RcppArmadillo")) install.packages("RcppArmadillo")
if (!require("MASS")) install.packages("MASS")
if (!require("dplyr")) install.packages("dplyr")

library(Rcpp)
library(MASS)
library(dplyr)

# ------------------------------------------------------------------------------
# Module 1: Embed the C++ Gibbs Sampler
# ------------------------------------------------------------------------------
# This C++ function implements a Bayesian Ridge Regression (BRR) Gibbs sampler
# optimized for multiple grouped variance components.
sourceCpp(code = '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List cpp_multikernel_gibbs(vec y, mat X, IntegerVector group_idx, int num_groups,
                           uvec missing_idx, uvec observed_idx,
                           int nIter, int burnIn) {
    int n = X.n_rows;
    int p = X.n_cols;
    
    // Initialize parameters
    vec beta = zeros<vec>(p);
    double mu = mean(y.elem(observed_idx));
    vec y_star = y - mu;
    
    // Variance components initialization
    vec var_u = ones<vec>(num_groups) * (var(y.elem(observed_idx)) / num_groups);
    double var_e = var(y.elem(observed_idx)) * 0.5;
    
    // Priors: Scaled Inverse Chi-Square
    double df = 5.0;
    double scale_e = var_e * (df - 2.0) / df;
    vec scale_u = var_u * (df - 2.0) / df;
    
    // Storage for Posterior Means
    vec post_beta = zeros<vec>(p);
    vec post_var_u = zeros<vec>(num_groups);
    double post_mu = 0.0;
    double post_var_e = 0.0;
    int samples = 0;
    
    // Precompute sum of squares for columns of X
    vec X_sq = sum(square(X), 0).t();
    uvec g_idx = as<uvec>(group_idx);
    
    for (int iter = 0; iter < nIter; ++iter) {
        
        // 1. Update Intercept (mu)
        double sum_e = sum(y_star.elem(observed_idx)) + mu * observed_idx.n_elem;
        double var_mu = var_e / observed_idx.n_elem;
        mu = randn<double>() * sqrt(var_mu) + (sum_e / observed_idx.n_elem);
        y_star = y - mu - X * beta;
        
        // 2. Update Marker Effects (Single-Site Gibbs)
        for (int j = 0; j < p; ++j) {
            int g = g_idx[j];
            double inv_var = (X_sq[j] / var_e) + (1.0 / var_u[g]);
            double mean_b = (dot(X.col(j), y_star) + X_sq[j] * beta[j]) / (var_e * inv_var);
            
            double old_beta = beta[j];
            beta[j] = randn<double>() * sqrt(1.0 / inv_var) + mean_b;
            y_star -= X.col(j) * (beta[j] - old_beta);
        }
        
        // 3. Update Variances (Gibbs via Inverse Gamma)
        double shape_e = (observed_idx.n_elem + df) / 2.0;
        double rate_e = (dot(y_star.elem(observed_idx), y_star.elem(observed_idx)) + df * scale_e) / 2.0;
        var_e = 1.0 / randg<double>(distr_param(shape_e, 1.0 / rate_e));
        
        for (int g = 0; g < num_groups; ++g) {
            uvec idx = find(g_idx == g);
            double shape_u = (idx.n_elem + df) / 2.0;
            double rate_u = (dot(beta.elem(idx), beta.elem(idx)) + df * scale_u[g]) / 2.0;
            var_u[g] = 1.0 / randg<double>(distr_param(shape_u, 1.0 / rate_u));
        }
        
        // 4. Impute Missing Phenotypes
        if (missing_idx.n_elem > 0) {
            vec y_hat_miss = mu + X.rows(missing_idx) * beta;
            y.elem(missing_idx) = y_hat_miss + randn<vec>(missing_idx.n_elem) * sqrt(var_e);
            
            // Re-sync y_star for the missing indices
            y_star.elem(missing_idx) = y.elem(missing_idx) - y_hat_miss;
        }
        
        // 5. Accumulate Posteriors
        if (iter >= burnIn) {
            post_beta += beta;
            post_var_u += var_u;
            post_mu += mu;
            post_var_e += var_e;
            samples++;
        }
    }
    
    return List::create(
        Named("mu") = post_mu / samples,
        Named("beta") = post_beta / samples,
        Named("var_u") = post_var_u / samples,
        Named("var_e") = post_var_e / samples
    );
}')


# ------------------------------------------------------------------------------
# Module 2: Data Simulation
# ------------------------------------------------------------------------------
simulate_example_met_data <- function(n_males = 15, n_females = 25, n_envs = 3, n_markers = 400) {
  set.seed(123)
  M_males <- matrix(sample(c(-1, 0, 1), n_males * n_markers, replace = TRUE), nrow = n_males)
  M_females <- matrix(sample(c(-1, 0, 1), n_females * n_markers, replace = TRUE), nrow = n_females)
  rownames(M_males) <- paste0("M", 1:n_males)
  rownames(M_females) <- paste0("F", 1:n_females)
  
  calc_G <- function(M) {
    W <- scale(M, center = TRUE, scale = FALSE)
    G <- tcrossprod(W) / (2 * sum(colMeans(M + 1)/2 * (1 - colMeans(M + 1)/2)))
    diag(G) <- diag(G) + 1e-4 
    return(G)
  }
  
  # Generate Unbalanced Data (30% missing)
  full_design <- expand.grid(Male = rownames(M_males), Female = rownames(M_females), Env = paste0("E", 1:n_envs))
  keep_idx <- sample(1:nrow(full_design), size = floor(nrow(full_design) * 0.7))
  unbalanced_data <- full_design[keep_idx, ]
  
  unbalanced_data$Pheno <- rnorm(nrow(unbalanced_data), 100, 15) # Dummy phenotypes
  
  return(list(observed_data = unbalanced_data, G_M = calc_G(M_males), G_F = calc_G(M_females)))
}

# ------------------------------------------------------------------------------
# Module 3: Eigen-Decomposition and Prediction Space Construction
# ------------------------------------------------------------------------------
build_eigen_space <- function(observed_data, G_M, G_F) {
  envs <- unique(observed_data$Env)
  full_data <- expand.grid(Male = rownames(G_M), Female = rownames(G_F), Env = envs) %>%
    left_join(observed_data, by = c("Male", "Female", "Env"))
  
  Z_M <- model.matrix(~ factor(Male, levels = rownames(G_M)) - 1, data = full_data)
  Z_F <- model.matrix(~ factor(Female, levels = rownames(G_F)) - 1, data = full_data)
  Z_E <- model.matrix(~ factor(Env, levels = envs) - 1, data = full_data)
  
  K_M <- Z_M %*% G_M %*% t(Z_M)
  K_F <- Z_F %*% G_F %*% t(Z_F)
  K_SCA <- K_M * K_F
  
  # Helper function to perform Eigen decomposition and return Z = U * sqrt(D)
  get_eigen_Z <- function(K) {
    eig <- eigen(K, symmetric = TRUE)
    # Keep only positive eigenvalues to avoid complex numbers
    pos_idx <- eig$values > 1e-8
    return(eig$vectors[, pos_idx] %*% diag(sqrt(eig$values[pos_idx])))
  }
  
  cat("Performing Eigen-decomposition on kernels...\n")
  Z_list <- list(
    Male_GCA = get_eigen_Z(K_M),
    Female_GCA = get_eigen_Z(K_F),
    SCA = get_eigen_Z(K_SCA)
  )
  
  # Concatenate matrices and build group indices for C++
  X_matrix <- do.call(cbind, Z_list)
  group_idx <- rep(0:(length(Z_list) - 1), times = sapply(Z_list, ncol))
  
  return(list(full_data = full_data, X = X_matrix, group_idx = group_idx, Z_dims = sapply(Z_list, ncol)))
}


# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

# 1. Setup
sim <- simulate_example_met_data()
prep <- build_eigen_space(sim$observed_data, sim$G_M, sim$G_F)

y_pheno <- prep$full_data$Pheno

# Create C++ 0-based indices for observed and missing data
missing_idx <- which(is.na(y_pheno)) - 1
observed_idx <- which(!is.na(y_pheno)) - 1

# Replace NAs with 0 before passing to C++
y_pheno[is.na(y_pheno)] <- 0

# 2. Run C++ Gibbs Sampler
cat("Running C++ Gibbs Sampler...\n")
fit_cpp <- cpp_multikernel_gibbs(
  y = y_pheno,
  X = prep$X,
  group_idx = prep$group_idx,
  num_groups = length(prep$Z_dims),
  missing_idx = missing_idx,
  observed_idx = observed_idx,
  nIter = 3000, 
  burnIn = 1000
)

# 3. Reconstruct GEBVs in R
cat("Reconstructing GEBVs...\n")

# Extract the beta segments corresponding to each kernel
idx_M <- 1:(prep$Z_dims[1])
idx_F <- (prep$Z_dims[1] + 1):(prep$Z_dims[1] + prep$Z_dims[2])
idx_SCA <- (prep$Z_dims[1] + prep$Z_dims[2] + 1):sum(prep$Z_dims)

# Re-project the betas back into the observation space (X * beta = u)
u_M <- prep$X[, idx_M] %*% fit_cpp$beta[idx_M]
u_F <- prep$X[, idx_F] %*% fit_cpp$beta[idx_F]
u_SCA <- prep$X[, idx_SCA] %*% fit_cpp$beta[idx_SCA]

prep$full_data$GEBV <- fit_cpp$mu + u_M + u_F + u_SCA

# 4. Summarize and Rank
top_crosses <- prep$full_data %>%
  select(Male, Female, GEBV) %>%
  distinct() %>%
  arrange(desc(GEBV)) %>%
  mutate(Rank = row_number())

cat("\n--- Top 10 Crosses (C++ Accelerated) ---\n")
print(head(top_crosses, 10))