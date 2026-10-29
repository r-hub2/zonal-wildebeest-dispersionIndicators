# R package `dispersionIndicators`
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!--  [![MetaboHUB Logo](https://forge.inrae.fr/pfem/gitlab-tools/-/raw/init/projects-templates/logos/metabohub_logo-20x20.png)![MetaboHUB title](https://img.shields.io/badge/MetaboHub-Software-0066cc?style=flat-square)](https://www.metabohub.fr)  -->
[![Latest Release](https://forge.inrae.fr/pfem/dev/libs/dispersion_indicators/-/badges/release.svg)](https://forge.inrae.fr/pfem/dev/libs/dispersion_indicators/-/releases)
[![coverage report](https://forge.inrae.fr/pfem/dev/libs/dispersion_indicators/badges/dev/coverage.svg)](https://forge.inrae.fr/pfem/dev/libs/dispersion_indicators/-/commits/dev)
[![pipeline status](https://forge.inrae.fr/pfem/dev/libs/dispersion_indicators/badges/dev/pipeline.svg)](https://forge.inrae.fr/pfem/dev/libs/dispersion_indicators/-/commits/dev)  

## Metadata

- authors: <etienne.jules@inrae.fr> 
- creation date: `2025-07-28`
- main usage: This repository contains the code of the R package dispersionIndicators.

## Description

This package provides function for analysing the dispersion of values from a 
batched and ordered dataset.

It can generate convex hull plots for each batch of a given dataset, along with
computing convex and ICM (Integrated Covariance Mahalanobis) dispersion
indicators.

The convex hull of a set of multidimensional points is the smallest convex
polygon that encloses all its points.

The convex indicators are based on convex hulls while  ICM dispersion indicators
are rather based on mahalanobis distance; All these provide insights into the
distribution of the data points and their dispersion among and between batches.

These indicators where developped in the context of metabolomics data during the 
thesis of Elfried Salanon [[1](#1)] that ended in 2025 (soon to be published). The convex
indicators are extensively described in [[2](#2)] while the ICM indicators are
detailed in [[3](#3)].

## Features

The main functions in this package are:
- `compute_icm_distances(...)`: Computes ICM 
    dispersion indicators.
- `convex_analysis_of_variables(...)`: Performs a comprehensive analysis by
    computing convex hull and all convex indicators for each
    variable (columns) in the dataset.
- `plot_all_convex_hulls(...)`: Generates convex hull plots
    for each batch in the dataset.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 
See deployment for notes on how to deploy the project on a live system.

### Prerequisites

<!-- NOTE: this section is required; list what things you need to install the software and how to install them -->

Before you begin, ensure you have met the following requirements:

- **Operating System**: Windows 10+, macOS 10.14+, or Linux
- **Programming Language**: [R] version 4.1 or higher
- **Other dependencies**: R packages ggplot2, stats, corpcor, utils

### Installing

You can install the package from CRAN using the following command:
```R
install.packages("dispersionIndicators")
```

### Running the tests


You can run the test to check eveything works fine as follows

 1. Open a R console
 2. (Assuming you have devtools installed) run:
  ```
  devtools::test()
  ```

This runs the unit test that are implemented in `dispersionIndicators/tests`

### Usage
Load the package:
```R
library(dispersionIndicators)
```
Prepare your data:
```R
data <- data.frame(
batch = rep(c("A","B","C"), each = 10),
injectionOrder = rep(1:30, times = 1),
metabolite1 = rnorm(30, mean = 100, sd = 10),
metabolite2 = rnorm(30, mean = 200, sd = 20)
)
```
Run the analysis:
```
result <- convex_analysis_of_variables(
  data = data,
  variable_columns = c("metabolite1", "metabolite2"),
  batch_col = "batch",
  sample_order_col = "injectionOrder",
  impute_if_needed = "median",
  mode = "global"
)
```
Generate and save the convex hull plots:
```
plot_all_convex_hulls(
  target_file_path = "../convex_hulls.pdf",
  convex_analysis_res = result,
  show_points = TRUE,
  mode = "global"
)
```
Compute ICM dispersion indicators:
```
result <- compute_icm_distances(
  data,
  batch_col = "batch",
  mode = "all",
  center_method_individual = "batch",
  center_method_inter = "mean"
)
```

## Changelog


All notable changes to this project will be documented in [CHANGELOG.md](NEWS.md).
For the versions available, see the [tags on this repository](https://forge.inrae.fr/pfem/dev/libs/dispersion_indicators/-/tags).

## Authors

- **Elfried Salanon** - *Theoretical work and idea, first draft of implementation* - MTH, INRAE, UNH, PFEM
- **Brice Mulot** - *First packaging and code cleaning pass* - INRAE, UNH, PFEM
- **Marie Lefebvre** - *Project management* - MTH, INRAE, UNH, PFEM
- **Etienne Jules** - *Final packaging, refactoring, project management, CRAN submission, conda forge submission, maintainer* - MTH, INRAE, UNH, PFEM

## License

**dispersionIndicators** is distributed under the ~~MIT License~~.
Please refer to [LICENSE](LICENSE) file for further details.

## Support &amp; External resources

- :book: **Documentation** - [Full documentation](https://CRAN.R-project.org/package=dispersionIndicators)
- :bug: **Bug Reports** - [GitLab Issues](https://forge.inrae.fr/pfem/dev/libs/dispersion_indicators/-/work_items)
- :email: **Email** - <etienne.jules@inrae.fr>

## Acknowledgments

<!-- NOTE: this section is required -->

- Thanks to Elfried Salanon, Brice Mulot, Marie Lefebvre and Etienne Jules for their valuable contributions

---

## References

<a name="1">[1]</a> **Salanon, Elfried M.B. (2025)**.
*Reproducibility and multi-source integration  in metabolomics for the
identification of common phenotypes.*

<a name="2">[2]</a> **Salanon, Elfried M.B. & al. (2024)**.
*An alternative for the robust assessment of the repeatability and
reproducibility of analytical measurements using bivariate dispersion.*
Chemometrics and Intelligent Laboratory Systems, 250, 105148.

<a names="3">[3]</a> **Salanon, Elfried M.B. & al. (2025)**.
*A benchmarking workflow for assessing the reliability of batch correction methods*