# HybridPred-MET: Genomic Prediction of Hybrid Performance in Multi-Environment Trials

An optimized R and C++ analytical pipeline for predicting hybrid crop performance across multi-environment trials (MET). This framework partitions additive (General Combining Ability, GCA), dominance (Specific Combining Ability, SCA), and genotype-by-environment (G×E) interaction variances using Bayesian Reproducing Kernel Hilbert Spaces (RKHS) regressions. 

## 🚀 Features

*   **Multi-Environment Architecture:** Natively models main genetic effects alongside $\text{GCA} \times E$ and $\text{SCA} \times E$ interactions using Hadamard products.
*   **Dual Solver Engine:**
    *   **BGLR Implementation:** Standard Gibbs sampling with built-in missing data imputation for T2 cross-validation schemes (predicting untested crosses from tested parents).
    *   **C++ Accelerated Solver:** A custom multi-threaded Bayesian Ridge Regression (BRR) engine built with `RcppArmadillo`. By utilizing Eigen-decomposition of the RKHS kernels ($K = U \Lambda U^T$), it reduces $O(N^3)$ matrix inversions to $O(N \times p)$ complexity, allowing seamless scaling to massive breeding datasets.
*   **Actionable Breeding Metrics:** Automated extraction of broad-adaptability Genomic Estimated Breeding Values (GEBVs) and calculation of plot-basis and entry-mean basis heritabilities ($h^2$ and $H^2$) for unbalanced designs.
*   **Diagnostic Visualizations:** Integrated functions to plot variance components, genetic effect distributions, and crossover interaction reaction norms.

## 🛠️ Prerequisites & Installation

The pipeline requires R and a C++ compiler. 

### Required R Packages
```R
install.packages(c("BGLR", "Rcpp", "RcppArmadillo", "MASS", "dplyr", "tidyr", "ggplot2", "gridExtra"))

```

## 📂 Pipeline Structure

The repository is structured into distinct modular stages:

1. **`01_Data_Simulation.R`**: Generates a synthetic, unbalanced MET dataset with predefined genetic architectures to validate the prediction pipeline prior to analyzing empirical data.
2. **`02_Kernel_Construction.R`**: Computes VanRaden $\mathbf{G}$ matrices and projects them into the full combinatorial factorial testing space (Male $\times$ Female $\times$ Env).
3. **`03_BGLR_Engine.R`**: Fits the model using the standard `BGLR` package, automatically imputing unobserved phenotypes.
4. **`04_Cpp_Accelerated_Engine.cpp` & `.R**`: The Eigen-decomposition pipeline and embedded C++ Gibbs sampler for rapid, large-scale dataset resolution.
5. **`05_Diagnostics_and_Heritability.R`**: Variance extraction, entry-mean heritability calculation handling unbalanced trial harmonic means, and `ggplot2` visualizations.

## 💻 Quick Start (C++ Accelerated Workflow)

```R
library(Rcpp)
library(dplyr)
sourceCpp("src/cpp_multikernel_gibbs.cpp")
source("R/pipeline_functions.R")

# 1. Load your data or run the simulation
sim <- simulate_example_met_data(n_males = 15, n_females = 25, n_envs = 3, prop_missing = 0.3)

# 2. Build the Eigen-decomposed prediction space
prep <- build_eigen_space(observed_data = sim$observed_data, G_M = sim$G_M, G_F = sim$G_F)

# 3. Format inputs for C++
y_pheno <- prep$full_data$Pheno
missing_idx <- which(is.na(y_pheno)) - 1
observed_idx <- which(!is.na(y_pheno)) - 1
y_pheno[is.na(y_pheno)] <- 0 # Initialize NAs

# 4. Run the RcppArmadillo solver
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

# 5. Extract GEBVs for selection
top_crosses <- extract_and_rank_gebvs(fit_cpp, prep)
print(head(top_crosses, 10))

```

## 📊 Output Interpretation

The framework outputs two primary targets for quantitative geneticists:

* **GEBV (Broad Adaptability):** The sum of the population mean ($\mu$) and the main genetic effects ($\text{GCA}_M + \text{GCA}_F + \text{SCA}$). Environmental deviations are excluded to isolate the core genetic potential of the cross across all target environments.
* **Entry-Mean Heritability ($H^2_{mean}$):** Adjusted for the unbalanced nature of the MET design using the harmonic mean of trial environments, providing a reliable metric for expected genetic gain.

## 📝 References

* Basnet, B. R., et al. (2019). "Modeling Genotype × Environment Interaction for Yield Prediction of Multi-Environment Trials." *The Plant Genome*, 12(2).
* Pérez, P., & de los Campos, G. (2014). Genome-wide regression and prediction with the BGLR statistical package. *Genetics*, 198(2), 483-495.

---

**Author / Maintainer:** Umesh R. Rosyara

**License:** MIT

```

```
