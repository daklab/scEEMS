# Minimal Training Demo (Chromosome 2 Only)

This is the public minimal model-training path using only chr2 train/test parquet files.

For full LOCO-style training across many chromosomes, download the full `model_training` resource from Synapse and run the standard command in `README.md`.

## Required files (chr2-only)

Place these files in your data bundle:

1. Train parquet (chr2):
- `training_data/Mic_mega_eQTL/training_data/train_NPR_10_PIP_0.1_0.01/annotated_data_Mic_mega_eQTL_chr2.parquet/`

2. Test parquet (chr2):
- `training_data/Mic_mega_eQTL/training_data/test_NPR_10_PIP_0.9_0.01/annotated_data_Mic_mega_eQTL_chr2.parquet/`

3. Supporting files:
- `<AUX_ROOT>/gnomad_MAF_chr2.tsv`
- `<AUX_ROOT>/columns_dict.pkl`
- `<AUX_ROOT>/41588_2024_1820_MOESM4_ESM.xlsx`


| File | Where it should live after download | Download URL                                                                                                                    |
|---|---|---------------------------------------------------------------------------------------------------------------------------------|
| `annotated_data_Mic_mega_eQTL_chr2.parquet` (train) | `training_data/Mic_mega_eQTL/training_data/train_NPR_10_PIP_0.1_0.01/` | `https://www.dropbox.com/scl/fo/6lvuvhj4d933yc30otvq4/ALgMp8ePqZEZMO_DKNOqiQI?rlkey=ch597t5je77q7a8nq6bws17g9&st=fymarlt2&dl=0` |
| `annotated_data_Mic_mega_eQTL_chr2.parquet` (test) | `training_data/Mic_mega_eQTL/training_data/test_NPR_10_PIP_0.9_0.01/` | `https://www.dropbox.com/scl/fo/q4c4aahxg5m7b2ltd1btf/AMidHs3CJFv611oryvFp06c?rlkey=rqzx18m46vgd8i4g2glu7tgpm&dl=0`             |
| `gnomad_MAF_chr2.tsv` | `<AUX_ROOT>/gnomad_MAF_chr2.tsv` | `https://www.dropbox.com/scl/fi/omlbn42b28moq6zhrw8lo/gnomad_MAF_chr2.tsv?rlkey=f4gqrzkl8apemldy84nmq6nga&dl=0`                 |
| `columns_dict.pkl` | `<AUX_ROOT>/columns_dict.pkl` | `https://www.dropbox.com/scl/fi/gwvqlejde8kmr1k8fjvvr/columns_dict.pkl?rlkey=yz8nln5pbbvt6e3l8fs8z6rfh&dl=0`                    |
| `41588_2024_1820_MOESM4_ESM.xlsx` | `<AUX_ROOT>/41588_2024_1820_MOESM4_ESM.xlsx` | `https://www.dropbox.com/scl/fi/z0rnxwz695wc1k956w9ho/41588_2024_1820_MOESM4_ESM.xlsx?rlkey=6on7rg3ulbxlndgy1clk3v8ll&dl=0`     |

## Important: parquet inputs are directories, not single files

`annotated_data_Mic_mega_eQTL_chr2.parquet` must be a directory containing parquet part files (for example `part.0.parquet`, `part.1.parquet`, etc.).

### 1. Create train/test parquet directories

```bash
mkdir -p <DEMO_ROOT>/training_data/Mic_mega_eQTL/training_data/train_NPR_10_PIP_0.1_0.01/annotated_data_Mic_mega_eQTL_chr2.parquet
mkdir -p <DEMO_ROOT>/training_data/Mic_mega_eQTL/training_data/test_NPR_10_PIP_0.9_0.01/annotated_data_Mic_mega_eQTL_chr2.parquet
```

### 2. Download all part files into each directory

- From the **train** Dropbox folder URL, download every file in that folder and place them in:
  - `<DEMO_ROOT>/training_data/Mic_mega_eQTL/training_data/train_NPR_10_PIP_0.1_0.01/annotated_data_Mic_mega_eQTL_chr2.parquet/`
- From the **test** Dropbox folder URL, download every file in that folder and place them in:
  - `<DEMO_ROOT>/training_data/Mic_mega_eQTL/training_data/test_NPR_10_PIP_0.9_0.01/annotated_data_Mic_mega_eQTL_chr2.parquet/`

Do not rename these directories or flatten their contents.

### 3. Quick check before training

```bash
find <DEMO_ROOT>/training_data/Mic_mega_eQTL/training_data/train_NPR_10_PIP_0.1_0.01/annotated_data_Mic_mega_eQTL_chr2.parquet -type f
find <DEMO_ROOT>/training_data/Mic_mega_eQTL/training_data/test_NPR_10_PIP_0.9_0.01/annotated_data_Mic_mega_eQTL_chr2.parquet -type f
```

## Minimal config snippet

```yaml
paths:
  data_dir: "<DEMO_ROOT>"
  gnomad_maf_dir: "<AUX_ROOT>"
  columns_dict_file: "<AUX_ROOT>/columns_dict.pkl"
  scratch_dir: "/tmp"
```

## Run minimal training (chr2-only)

Use the wrapper (recommended):

```bash
cd 5_model_training
bash run_minimal_training.sh \
  Mic_mega_eQTL \
  2 \
  "<AUX_ROOT>/41588_2024_1820_MOESM4_ESM.xlsx"
```

Equivalent direct command:

```bash
cd 5_model_training
python train_model.py Mic_mega_eQTL 2 \
  --gene_lof_file "<AUX_ROOT>/41588_2024_1820_MOESM4_ESM.xlsx" \
  --yaml_path data_params.yaml \
  --single_chromosome_demo
```

## Important note

- `--single_chromosome_demo` is for public minimal demonstration only.
- It is not the full multi-chromosome LOCO training setup used for the main production workflow.
