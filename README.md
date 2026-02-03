# DNA Methylation Mixture Deconvolution Scripts

This repository contains the R scripts used for the analysis of two-person DNA methylation mixtures in our forensic epigenetics study:

**Poggiali et al., Forensic Science International: Genetics (2026)**  
DOI: https://doi.org/10.1016/j.fsigen.2026.103438

These scripts were developed to:
- estimate contributor-specific methylation profiles from mixed samples
- support mixture deconvolution experiments
- benchmark deconvoluted profiles against single-source references
- reproduce key analyses reported in the manuscript

The repository is intentionally minimal and contains only the functions required to reproduce the published results.

---

## 📦 Scripts:
- **EPICv2.0_analysis_deconvolution_and_age_pred_DNAm_mixture.R** → Analysis of the DNA mixtures generated in the laboratory  
- **Generation_and_deconvolution_in_silico_generated_mixture.R** → Analysis of the in-silico DNA mixtures  

---

## ⭐ UnMixMe (UnMix DNA Methylation profiles)

The following function performs deconvolution of the offender DNA methylation profile from a two-person mixture (UnMixMe approach):

```r
mixture_deconvolution <- function(beta_mixture, beta_victim,
                                  proportion_victim = 1,
                                  proportion_suspect = 1){

  beta_offender <- (((proportion_victim + proportion_suspect) * beta_mixture) -
                    (proportion_victim * beta_victim)) / proportion_suspect

  return(beta_suspect)
}
