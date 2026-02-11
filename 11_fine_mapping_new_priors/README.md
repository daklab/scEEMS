# Step 11: Fine-Mapping with scEEMS Priors

Re-run SuSiE statistical fine-mapping using scEEMS predictions as informed priors.

## Overview

This step demonstrates how scEEMS predictions can improve statistical fine-mapping by providing biologically-informed prior probabilities for each variant being causal. The approach:

1. Calibrates scEEMS prediction probabilities to a genomic base rate (1/1000)
2. Runs SuSiE fine-mapping with these calibrated priors
3. Compares results to standard SuSiE with uniform priors

## Scripts

| Script | Description |
|--------|-------------|
| `download_data.py` | Download reference panel from Synapse (requires token in config.yaml) |
| `template_predicted.R` | SuSiE fine-mapping with scEEMS priors |
| `template_uniform.R` | Baseline SuSiE with uniform priors |
| `compare_pips.py` | Compare PIPs between uniform and predicted priors |
| `run_pipeline.sh` | Run full fine-mapping comparison |

## Prior Calibration

scEEMS prediction probabilities are converted to fine-mapping priors using logit adjustment:
```
logit_prior = log(pred_prob / (1 - pred_prob))
logit_modified = logit_prior + log(base_rate / (1 - base_rate) / (0.5 / 0.5))
calibrated_prior = 1 / (1 + exp(-logit_modified))
```

Where `base_rate = 1/1000` is the expected proportion of causal variants.

## Inputs

- scEEMS predictions from Step 6
- Genotype, phenotype, and covariate files
- Reference panel (Synapse: syn53171227)

## Outputs

- `{fine_mapping_dir}/{cell_type}/fine_mapping_predicted/` - RDS files with scEEMS priors
- `{fine_mapping_dir}/{cell_type}/fine_mapping_uniform/` - RDS files with uniform priors
- Comparison TSV files
