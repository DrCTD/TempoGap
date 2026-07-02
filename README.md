# TempoGap

TempoGap is an R package for estimating the timing of immune recovery from transcriptomic or proteomic profiles. The algorithm trains immune recovery clocks to predict biological time after injury or critical illness and then calculates the temporal difference between predicted biological recovery time and actual clinical time.

The resulting **TempoGap score** can be used to classify samples into biologically **Leading**, **Intermediate**, or **Lagging** recovery states.

---

## Overview

TempoGap was developed to quantify deviations in temporal immune recovery among patients with trauma or critical illness.

The general workflow is:

```text
Expression matrix + clinical/sample metadata
        ↓
Train immune recovery clock
        ↓
Predict biological recovery time
        ↓
Calculate TempoGap = predicted time - actual time
        ↓
Classify samples into Leading / Intermediate / Lagging recovery states
```

A positive TempoGap indicates that the sample appears biologically more advanced in recovery than expected for its actual clinical time.  
A negative TempoGap indicates that the sample appears biologically delayed relative to its actual clinical time.

---

## Installation

TempoGap can be installed directly from GitHub:

```r
install.packages("devtools")
devtools::install_github("DrCTD/TempoGap")
```

Load the package:

```r
library(TempoGap)
```

---

## Tutorial

A step-by-step HTML tutorial is provided to help users train their own immune recovery clocks and calculate TempoGap scores.

Tutorial link:

[TempoGap vignette](https://github.com/DrCTD/TempoGap/blob/main/vignette.html)

The tutorial introduces how to:

```text
1. Load expression and metadata files
2. Train an immune recovery clock
3. Predict biological recovery time
4. Calculate TempoGap scores
5. Classify samples into recovery states
6. Visualize TempoGap results
```

---

## Example data

Example data are available through the GitHub Release page:

[Download TempoGap example data](https://github.com/DrCTD/TempoGap/releases/tag/TempoGap)

After downloading the example data, unzip the file and follow the HTML tutorial.

The example data are provided for demonstration purposes and can be used to test the TempoGap workflow.

---

## Input format

### Expression matrix

The expression matrix should contain genes or proteins as rows and samples as columns.

```text
Gene/Protein    Sample1    Sample2    Sample3
ANXA6           7.21       6.98       7.55
IL4R            5.31       5.89       6.12
IL1R2           8.04       7.76       8.51
```

### Metadata

The metadata table should contain samples as rows and clinical or sample-level variables as columns.

```text
sample_id    hours_since_injury    recovery_class
Sample1      12                    UCR
Sample2      72                    IR
Sample3      144                   CR
```

The sample names in the metadata must match the column names of the expression matrix.

Before running TempoGap, we recommend checking sample matching:

```r
common_samples <- intersect(colnames(expr), rownames(meta))

expr <- expr[, common_samples]
meta <- meta[common_samples, ]

all(colnames(expr) == rownames(meta))
```

---

## Quick start

```r
library(TempoGap)

# Load your expression matrix and metadata
expr <- read.csv(
  "example_expression_matrix.csv",
  row.names = 1,
  check.names = FALSE
)

meta <- read.csv(
  "example_metadata.csv",
  row.names = 1,
  check.names = FALSE
)

# Train a TempoGap model
model <- train_temporal_model(
  expr_mat = as.matrix(expr),
  metadata = meta,
  time_col = "hours_since_injury"
)

# Predict TempoGap scores
tempogap_result <- predict_temporal_gap(
  model = model,
  expr_mat = as.matrix(expr),
  metadata = meta,
  time_col = "hours_since_injury"
)

head(tempogap_result)
```

---

## Output

The TempoGap output includes:

```text
sample_id
actual_time
predicted_time
TempoGap
TempoGap_z
TempoGap_group
```

TempoGap is calculated as:

```text
TempoGap = predicted biological time - actual clinical time
```

A typical classification strategy is:

```text
TempoGap_z > 1      = Leading
TempoGap_z < -1     = Lagging
Otherwise           = Intermediate
```

Users may adjust the threshold according to their study design and analysis goals.

---

## Main functions

| Function                           | Description                                                  |
| ---------------------------------- | ------------------------------------------------------------ |
| `train_temporal_model()`           | Trains an immune recovery clock using expression data and actual clinical time |
| `predict_temporal_gap()`           | Predicts biological time and calculates TempoGap scores      |
| `plot_bootstrap_coefficients()`    | Trains bootstrapped temporal models                          |
| `plot_temporal_gap_distribution()` | Visualizes TempoGap score distribution                       |
| `cox_analysis_tempo_gap()`         | Tests associations between TempoGap and time-to-event outcomes |
| `auc_analysis_tempo_gap()`         | Evaluates classification performance using ROC/AUC analysis  |

---

## Data availability

The example data provided in this repository are intended for demonstration and tutorial use.

Full transcriptomic, proteomic, and clinical datasets associated with the study are available through the repositories described in the manuscript Data Availability section.

---

## Code availability

The TempoGap algorithm is available as an R package in this GitHub repository. The repository includes source code, documentation, an HTML tutorial, and example data to reproduce the main TempoGap workflow.

---

## Contact

For questions, please contact:

```text
Teding Chang
changt5@upmc.edu
Department of Surgery
University of Pittsburgh
```

---

## License

This package is released under the license specified in the `LICENSE` file.