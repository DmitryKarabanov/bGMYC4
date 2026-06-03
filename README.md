# bGMYC4
## *Reid's Bayesian species delimitation, made reproducible and accessible*

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**Bayesian General Mixed Yule-Coalescent (GMYC) model for species delimitation in R ≥ 4.5**

---

## 🌍 Overview 

bGMYC4 implements a Bayesian version of the General Mixed Yule-Coalescent model for species delimitation from single-locus data. It uses Markov Chain Monte Carlo (MCMC) simulation to evaluate the posterior distribution of species boundaries, allowing the use of multiple phylogenetic trees to account for phylogenetic uncertainty. The package includes interactive diagnostics, parameter tuning, and visualization tools. Full documentation and interactive vignette available in both English and Russian. Full simple interactiv scrip available in inst\scripts

---

## ✨ Key Features 

| Feature | Description |
|---------|-------------|
| 🌲 Single & Multi-tree analysis | Run GMYC on one tree or pool uncertainty across multiple BEAST2 trees |
| 📊 Interactive diagnostics | Real-time MCMC checks: trace plots, acceptance rates, ESS, Gelman-Rubin R̂, auto-tuning hints |
| ⚙️ Parameter tuning | Safe input handlers with biological hints for priors, scales, and thresholds |
| 🌡️ Heatmap visualization | Conspecificity probability matrices ordered by tree topology |
| 🌐 Bilingual workflow | Full EN/RU documentation and interactive vignette |
| 🧬 BEAST2 ready | Automatic cleaning of Nexus annotations and ultrametricity enforcement |
| ⚡ Optimized performance | Parallel `bgmyc.multiphylo()`, pure-R post-processing, vectorized matrices |
| ✅ Full CRAN compliance | `0 errors, 0 warnings, 0 notes`|
| ⚡ High-performance prallelization | Stable parallel backend for Windows |
| 🔍 Detailed guide | Interactive vignette with dynamic convergence diagnostics (need knitr, rmarkdown, pandoc). Full simple interactiv scrip available in `inst\scripts` |
| 🧹 Pure-R post-processing | No external `sort`/`uniq` dependencies |

---

## 📦 Installation 

### 🔹 From GitHub (recommended)
```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("DmitryKarabanov/bGMYC4", build_vignettes = TRUE)
```

### 🔹 From source archive
```r
install.packages("path/to/bGMYC4_4.1.0.tar.gz", repos = NULL, type = "source")
# To build vignettes locally:
# install.packages("path/to/bGMYC4.tar.gz", repos = NULL, type = "source", build_vignettes = TRUE)
```

### 🔹 System requirements
- **R ≥ 4.5.1** (tested on 4.5.1–4.5.3)
- **OS:** Windows 11 x64, macOS, Linux
- **Dependencies:** `ape`, `future`, `future.apply`, `parallel`, `knitr`, `rmarkdown`
- **Recommended for Diagnostics:** `mcmcse` (for ESS & Gelman-Rubin $\hat{R}$ convergence checks)

---

## 🚀 Quick Start 

### 1. Load the package
```r
library(bGMYC4)
```

### 2. Open the interactive tutorial
```r
browseVignettes("bGMYC4")
# Or run directly:
vignette("bGMYC4-interactive", package = "bGMYC4")
```

### 3. Run a minimal example
```r
# Load built-in test data
data(east10)

# Run single-tree analysis (short MCMC for testing)
result <- bgmyc.singlephy(
  phylo = east10[[1]],
  mcmc = 1000, burnin = 200, thinning = 10,
  py1 = 0, py2 = 1.5, pc1 = 0, pc2 = 2,
  t1 = 2, t2 = 35,
  scale = c(20, 10, 5),
  start = c(1, 0.5, 10)
)

# Check convergence
plot(result)

# Build conspecificity matrix
prob_mat <- spec.probmat(result)
plot(prob_mat, east10[[1]])

# Extract species clusters (your threshold)
out <- bgmyc.point(prob_mat, ppcutoff = 0.05)
cat("Delimited clusters:", length(out), "\n")
```

---

## 📚 Documentation 

| Resource | Command | Description |
|----------|---------|-------------|
| Package help | `?bGMYC4` | Function reference |
| Interactive tutorial | `browseVignettes("bGMYC4")` | Step-by-step workflow with diagnostics |
| Function details | `?bgmyc.multiphylo`, `?spec.probmat` | Parameter explanations |
| Built-in data | `?east10` | Simulated dataset for testing |
| Performance tips | `vignette("bGMYC4-interactive", package = "bGMYC4")` | Section "Performance & Export" |

---

## 🔧 Workflow Overview 

1. **Data selection**  
   - Built-in `east10` (simulated)  
   - Custom BEAST2 files (`.tree`, `.trees`)

2. **Tree preprocessing**  
   - Strip BEAST2 annotations `[&...]`  
   - Rooting check and outgroup removal  
   - Ultrametricity enforcement (`fix_ultrametric()`)

3. **Interactive diagnostics**  
   - Run `bgmyc.singlephy()` on one tree  
   - Check acceptance rates and ESS  
   - Evaluate trace plots (`plot(result)`)  
   - Tune parameters with biological hints

4. **Multi-tree analysis**  
   - Run `bgmyc.multiphylo()` on sampled trees  
   - Pool phylogenetic uncertainty (parallel by default)

5. **Visualization and delimitation**  
   - Build conspecificity heatmap (`spec.probmat()`)  
   - Extract clusters at your `ppcutoff` (e.g., `0.05`)

---

## ⚙️ Parameter Reference 

| Parameter | Role | Recommended Range | Notes |
|-----------|------|------------------|-------|
| `mcmc` | Chain length | `10k` (test), `50k+` (final) | Longer = better posterior resolution |
| `burnin` | Burn-in | `20–30%` of `mcmc` | Discards non-stationary start |
| `thinning` | Sampling interval | `10–50` | Reduces autocorrelation & RAM usage |
| `py1`, `py2` | Yule rate priors | `0`, `0.5–1.5` | Model: λ ∝ n^py; >1.5 blurs Yule/Coalescent boundary |
| `pc1`, `pc2` | Coalescent priors | `0`, `1.0–2.0` | Models Ne change; <1 → decline, >1 → growth |
| `t1`, `t2` | Threshold prior (species count) | `2`, `min(y, ntips-1)` | Must be `< ntips`; auto-capped to prevent crashes |
| `scale` | MCMC proposal widths | `c(20–30, 10–15, 3–7)` | Tune via acceptance rates; higher = more conservative |
| `ppcutoff` | Species lumping threshold | `0.05` , `0.01` | Low = captures high intraspecific variation / ILS |

---

## ⚡ Performance & Compatibility 

The optimized `bGMYC4 v4.1.0` is **fully compatible with R ≥ 4.5** and includes:

| Feature / Возможность | Description |
|----------------------|----------------------|
| 🔁 **Automatic parallelization** | `bgmyc.multiphylo()` runs in parallel via `future.apply` (uses all CPU cores by default) |
| 🧹 **Pure-R post-processing** | `bgmyc.spec()` and `spec.probmat()` use native R functions — no external `sort`/`uniq` dependencies |
| 🚀 **Vectorized operations** | Faster heatmap construction and matrix indexing |
| 🔐 **Parallel-safe RNG** | `future.seed = TRUE` ensures reproducible multi-tree analyses |
| 🖥️ **R ≥ 4.5 compatibility** | Full support for modern R features, byte-compilation (`compiler::enableJIT(3)`), and optimized memory management |

### 📈 Statistical Rigor & Diagnostics 
- **Effective Sample Size (ESS)**: Automatically calculated per parameter when `mcmcse` is installed. Target: `ESS > 200`.
- **Gelman-Rubin $\hat{R}$**: Post-analysis check comparing MCMC chains across sampled trees. Automatically computed if `mcmcse` is available. Target: `R̂ < 1.05`.

### Limit CPU cores (optional) 
```r
# Before running bgmyc.multiphylo():
future::plan(future::multisession, workers = 8)  # Use 8 cores instead of all
```

### Export results to CSV 
```r
# After bgmyc.spec():
write.csv(spec_out$specprobs, 
          file = "bGMYC_clusters.csv", 
          row.names = FALSE, 
          fileEncoding = "UTF-8")
```

### Validate model fit 
```r
# Compare Yule vs Coalescent rates:
rates <- checkrates(result_multi)
plot(rates)  # Opens 4 diagnostic plots
```

### Speed tip 
```r
# Add to your ~/.Rprofile for automatic JIT acceleration:
compiler::enableJIT(3)  # +2–3% speedup on all MCMC loops, no code changes
```

---

## ❓ Troubleshooting 

| Issue | Solution |
|-------|----------|
| `Your input tree is not ultrametric` | Use `fix_ultrametric()` or check BEAST2 export settings |
| `start[3] out of bounds` | Ensure `t1 < start[3] < t2`; script auto-clamps values |
| `Error in read.nexus` | Verify files are Nexus format; use `read.tree()` for Newick |
| Non-ASCII warning in check | All R code now uses ASCII-only; comments may contain UTF-8 |
| `Error: package 'ape' is not available` | Run `install.packages("ape")` first |
| `future::evalFuture() failed` | Ensure package is installed |
| `Non-ASCII characters in Rd file` | All R code is ASCII; comments may contain UTF-8 — this is CRAN-compliant |
| `ESS: package 'mcmcse' not installed` | Run `install.packages("mcmcse")`. The package works without it, but full convergence stats require it. |

---

## 📄 Citation 

**If you use bGMYC4 in your research, please cite:**

```bibtex
@Manual{bGMYC4,
  title = {bGMYC4: Bayesian General Mixed Yule-Coalescent Model for Species Delimitation},
  author = {Karabanov, Dmitry},
  year = {2026},
  note = {R package version 4.1.0. Developed with assistance from Qwen3.6 AI assistant.},
  url = {https://github.com/DmitryKarabanov/bGMYC4},
 }
```

**Please, cite the foundational works:**
- Pons et al. (2006) *Syst. Biol.* 55:595 — GMYC model
- Reid & Carstens (2012) *Mol. Ecol. Res.* 12:446 — Bayesian implementation

---

## 🤝 Contributing 

Contributions are welcome! 

---

## 📜 License 

bGMYC4 is free software: you can redistribute it and/or modify it under the terms of the **GNU General Public License v3.0**.  
See [LICENSE](LICENSE) for details.

---

## 📬 Contact

**Dmitry Karabanov**  
📧 dk[at]ibiw.ru  
🏛 Institute for Biology of Inland Waters, RAS  
🔗 [ORCID](https://orcid.org/0000-0001-6008-7441)  

---

### 🤖 AI Assistance
ℹ️ This package was developed with assistance from Qwen3.6, an AI language model by Alibaba Cloud. The AI provided code optimization suggestions, debugging support, and documentation drafting. All scientific decisions, parameter choices, and biological interpretations were made by the author.    
[![Powered by Qwen](https://img.shields.io/badge/Powered%20by-Qwen-6366f1?style=flat-square&logo=data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMjQiIGhlaWdodD0iMjQiIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0ibm9uZSIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KPHBhdGggZD0iTTEyIDJMMTQuNTkgOC4yNkwyMSA5LjI3TDE2LjE4IDEzLjk3TDE3LjM0IDIwLjczTDEyIDE3LjU0TDYuNjYgMjAuNzNMNy44MiAxMy45N0wzIDkuMjdMOS40MSA4LjI2TDEyIDJaIiBmaWxsPSIjZmZmIi8+Cjwvc3ZnPg==)](https://qwen.ai/)

---
> *bGMYC4 – Bayesian species delimitation, made simple, reproducible, and biologist‑friendly* 🌿🔬
