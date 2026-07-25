# ==============================================================================
# MODULAR HYBRID PREDICTION PIPELINE (BGLR)
# General Combining Ability (GCA), Specific Combining Ability (SCA), and GxE
# ==============================================================================

if (!require("BGLR")) install.packages("BGLR")
if (!require("MASS")) install.packages("MASS")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("gridExtra")) install.packages("gridExtra")

library(BGLR)
library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# ------------------------------------------------------------------------------
# Module 1: Example Data Simulation
# Note: This generates example data for a user to work on if nothing else is available.
# ------------------------------------------------------------------------------
simulate_example_met_data <- function(n_males = 15, n_females = 25, n_envs = 3, 
                                      n_markers = 400, prop_missing = 0.3, seed = 123) {
  set.seed(seed)
  cat("Generating example dataset...\n")
  
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
  
  G_M <- calc_G(M_males)
  G_F <- calc_G(M_females)
  
  full_design <- expand.grid(Male = rownames(M_males), Female = rownames(M_females), Env = paste0("E", 1:n_envs))
  
  # Introduce missingness for unbalanced design
  n_total <- nrow(full_design)
  keep_idx <- sample(1:n_total, size = floor(n_total * (1 - prop_missing)))
  unbalanced_data <- full_design[keep_idx, ]
  rownames(unbalanced_data) <- NULL
  
  # Simulate effects
  n_obs <- nrow(unbalanced_data)
  Z_M <- model.matrix(~ factor(Male, levels = rownames(G_M)) - 1, data = unbalanced_data)
  Z_F <- model.matrix(~ factor(Female, levels = rownames(G_F)) - 1, data = unbalanced_data)
  Z_E <- model.matrix(~ factor(Env, levels = paste0("E", 1:n_envs)) - 1, data = unbalanced_data)
  
  K_M <- Z_M %*% G_M %*% t(Z_M)
  K_F <- Z_F %*% G_F %*% t(Z_F)
  K_SCA <- K_M * K_F                  
  K_E <- Z_E %*% t(Z_E)               
  
  mu <- 100
  env_effects <- mvrnorm(1, mu = rep(0, n_obs), Sigma = 10.0 * K_E)
  true_gca_m <- mvrnorm(1, mu = rep(0, n_obs), Sigma = 2.0 * K_M)
  true_gca_f <- mvrnorm(1, mu = rep(0, n_obs), Sigma = 2.5 * K_F)
  true_sca   <- mvrnorm(1, mu = rep(0, n_obs), Sigma = 1.5 * K_SCA)
  
  true_gxe_m <- mvrnorm(1, mu = rep(0, n_obs), Sigma = 1.0 * (K_M * K_E))
  true_gxe_f <- mvrnorm(1, mu = rep(0, n_obs), Sigma = 1.0 * (K_F * K_E))
  true_gxe_s <- mvrnorm(1, mu = rep(0, n_obs), Sigma = 0.5 * (K_SCA * K_E))
  
  residuals  <- rnorm(n_obs, 0, sqrt(3.0))
  
  unbalanced_data$Pheno <- mu + env_effects + true_gca_m + true_gca_f + 
    true_sca + true_gxe_m + true_gxe_f + true_gxe_s + residuals
  
  return(list(observed_data = unbalanced_data, G_M = G_M, G_F = G_F))
}

# ------------------------------------------------------------------------------
# Module 2: Kernel Construction (Full Factorial Space)
# ------------------------------------------------------------------------------
build_predictive_space <- function(observed_data, G_M, G_F) {
  cat("Building covariance kernels for the complete factorial space...\n")
  
  envs <- unique(observed_data$Env)
  full_design <- expand.grid(Male = rownames(G_M), Female = rownames(G_F), Env = envs)
  
  full_data <- full_design %>%
    left_join(observed_data %>% select(Male, Female, Env, Pheno), by = c("Male", "Female", "Env"))
  
  Z_M <- model.matrix(~ factor(Male, levels = rownames(G_M)) - 1, data = full_data)
  Z_F <- model.matrix(~ factor(Female, levels = rownames(G_F)) - 1, data = full_data)
  Z_E <- model.matrix(~ factor(Env, levels = envs) - 1, data = full_data)
  
  K_M <- Z_M %*% G_M %*% t(Z_M)
  K_F <- Z_F %*% G_F %*% t(Z_F)
  K_SCA <- K_M * K_F
  K_E <- Z_E %*% t(Z_E)
  
  kernels <- list(
    Env = list(X = Z_E, model = "BRR"), 
    Male_GCA = list(K = K_M, model = "RKHS"),
    Female_GCA = list(K = K_F, model = "RKHS"),
    SCA = list(K = K_SCA, model = "RKHS"),
    Male_Env = list(K = K_M * K_E, model = "RKHS"),
    Female_Env = list(K = K_F * K_E, model = "RKHS"),
    SCA_Env = list(K = K_SCA * K_E, model = "RKHS")
  )
  
  return(list(full_data = full_data, ETA = kernels))
}

# ------------------------------------------------------------------------------
# Module 3: Model Fitting & GEBV Extraction
# ------------------------------------------------------------------------------
fit_and_extract_gebvs <- function(full_data, ETA, nIter = 3000, burnIn = 1000) {
  cat(sprintf("\nFitting BGLR model (%d iterations)...\n", nIter))
  
  fit <- BGLR(y = full_data$Pheno, ETA = ETA, nIter = nIter, burnIn = burnIn, verbose = FALSE)
  
  full_data$GEBV <- fit$mu + fit$ETA$Male_GCA$u + fit$ETA$Female_GCA$u + fit$ETA$SCA$u
  
  hybrid_ranks <- full_data %>%
    select(Male, Female, GEBV) %>%
    distinct() %>%
    arrange(desc(GEBV)) %>%
    mutate(Rank = row_number()) %>%
    relocate(Rank)
  
  return(list(fit = fit, ranks = hybrid_ranks, full_data = full_data))
}

# ------------------------------------------------------------------------------
# Module 4: Entry-Mean Heritability Calculation
# ------------------------------------------------------------------------------
calc_heritability <- function(fit, observed_data) {
  # Harmonic mean of environments per cross
  env_counts <- observed_data %>% group_by(Male, Female) %>% summarize(n_env = n_distinct(Env), .groups = "drop")
  e_harmonic <- nrow(env_counts) / sum(1 / env_counts$n_env)
  r_harmonic <- 1 # Assuming 1 rep per environment
  
  V_A <- fit$ETA$Male_GCA$varU + fit$ETA$Female_GCA$varU
  V_G <- V_A + fit$ETA$SCA$varU
  V_GxE <- fit$ETA$Male_Env$varU + fit$ETA$Female_Env$varU + fit$ETA$SCA_Env$varU
  V_e <- fit$varE
  
  V_P_mean <- V_G + (V_GxE / e_harmonic) + (V_e / (e_harmonic * r_harmonic))
  
  cat(sprintf("\nEntry-Mean Heritabilities:\n Narrow-sense (h2): %.3f\n Broad-sense (H2):  %.3f\n", 
              V_A / V_P_mean, V_G / V_P_mean))
}

# ------------------------------------------------------------------------------
# Module 5: Visualization
# ------------------------------------------------------------------------------
plot_model_diagnostics <- function(fit, full_data) {
  # 1. Variance Components
  comp_names <- names(fit$ETA)
  var_list <- sapply(comp_names, function(x) ifelse(fit$ETA[[x]]$model == "BRR", fit$ETA[[x]]$varB, fit$ETA[[x]]$varU))
  var_list["Residual"] <- fit$varE
  
  p1 <- ggplot(data.frame(Comp = names(var_list), Var = var_list), aes(x = reorder(Comp, Var), y = Var)) +
    geom_bar(stat = "identity", fill = "#2c3e50") + coord_flip() + theme_minimal() +
    labs(title = "Variance Components", x = "", y = "Variance")
  
  # 2. Main Effects Distribution
  df_main <- data.frame(Male_GCA = fit$ETA$Male_GCA$u, Female_GCA = fit$ETA$Female_GCA$u, SCA = fit$ETA$SCA$u) %>% 
    pivot_longer(everything(), names_to = "Effect", values_to = "Value")
  
  p2 <- ggplot(df_main, aes(x = Value, fill = Effect)) +
    geom_density(alpha = 0.6) + scale_fill_brewer(palette = "Set1") + theme_minimal() +
    labs(title = "Genetic Effects Distribution", x = "Effect Size", y = "Density")
  
  # 3. GxE Reaction Norms (Top 5 Males)
  full_data$GxE_M <- fit$ETA$Male_Env$u
  full_data$GCA_M <- fit$ETA$Male_GCA$u
  
  top_males <- full_data %>% group_by(Male) %>% summarize(GCA = mean(GCA_M)) %>% top_n(5, GCA) %>% pull(Male)
  gxe_summary <- full_data %>% filter(Male %in% top_males) %>% group_by(Male, Env) %>% summarize(GxE = mean(GxE_M), .groups="drop")
  
  p3 <- ggplot(gxe_summary, aes(x = Env, y = GxE, group = Male, color = Male)) +
    geom_line(linewidth = 1) + geom_point(size = 3) + geom_hline(yintercept = 0, linetype = "dashed") +
    theme_minimal() + labs(title = "GxE Reaction Norms (Top 5 Males)", x = "Environment", y = "Deviation")
  
  grid.arrange(p1, p2, p3, ncol = 2, layout_matrix = rbind(c(1, 2), c(3, 3)))
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

# 1. Generate example data
sim_data <- simulate_example_met_data(n_males = 15, n_females = 25, n_envs = 3, prop_missing = 0.3)

# 2. Build kernels for the full factorial prediction space
prep <- build_predictive_space(sim_data$observed_data, sim_data$G_M, sim_data$G_F)

# 3. Fit model and extract GEBVs
results <- fit_and_extract_gebvs(prep$full_data, prep$ETA, nIter = 3000, burnIn = 1000)

# 4. Diagnostics & Reporting
calc_heritability(results$fit, sim_data$observed_data)
cat("\nTop 10 Crosses (Including Imputed):\n")
print(head(results$ranks, 10))

# 5. Plot (Uncomment to render)
# plot_model_diagnostics(results$fit, results$full_data)

# ==============================================================================
# CLEAN BGLR SPACE
# ==============================================================================
clean_bglr_traces <- function(prefix = "") {
  # Define the patterns BGLR uses (adjust prefix if you used 'saveAt' with a custom name)
  bglr_patterns <- paste0(
    "^", prefix, "ETA_.*\\.dat$|",  # ETA trace files (varU, varB, etc.)
    "^", prefix, "varE\\.dat$|",    # Residual variance
    "^", prefix, "mu\\.dat$"        # Intercept
  )
  
  # Find matching files in the current working directory
  files_to_remove <- list.files(pattern = bglr_patterns, full.names = TRUE)
  
  if (length(files_to_remove) > 0) {
    file.remove(files_to_remove)
    cat(sprintf("Successfully removed %d BGLR trace files.\n", length(files_to_remove)))
  } else {
    cat("No BGLR trace files found in the current directory.\n")
  }
}

# Execute the cleanup
clean_bglr_traces()
