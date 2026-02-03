# DNA Methylation Mixture Deconvolution Scripts

This repository contains lightweight R scripts used for the analysis of two-person DNA methylation mixtures in our forensic epigenetics study.

The code implements the core computational steps described in:

**Poggiali et al., Forensic Science International: Genetics (2026)**  
DOI: https://doi.org/10.1016/j.fsigen.2026.103438

These scripts were developed to:
- estimate contributor-specific methylation profiles from mixed samples
- support mixture deconvolution experiments
- benchmark deconvoluted profiles against single-source references
- reproduce key analyses reported in the manuscript

The repository is intentionally minimal and contains only the functions required to reproduce the published results.

---

## 📦 Repository structure

.
├── scripts/
│ ├── EPICv2.0_analysis_deconvolution_and_age_pred_DNAm_mixture.R
│ ├── Generation_and_deconvolution_in_silico_generated_mixture.R
└── README.md

- **EPICv2.0_analysis_deconvolution_and_age_pred_DNAm_mixture.R** → Analysis of the DNA mixtures generated in the laboratory  
- **Generation_and_deconvolution_in_silico_generated_mixture.R** → Analysis of the in-silico DNA mixtures  

---

## ⭐ Core function (used throughout the paper)

The following function performs deconvolution of the offender DNA methylation profile from a two-person mixture:

```r
mixture_deconvolution <- function(beta_mixture, beta_victim,
                                  proportion_victim = 1,
                                  proportion_offender = 1){

  beta_offender <- (((proportion_victim + proportion_offender) * beta_mixture) -
                    (proportion_victim * beta_victim)) / proportion_offender

  return(beta_offender)
}
