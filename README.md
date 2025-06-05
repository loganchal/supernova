# supernova

This repository contains the data analysis for the **Red Tide Supernova**
directed studies project.

## Installing Required Packages

Before running any of the analyses make sure the following R packages are
installed:

```R
install.packages(c("tidyverse", "ggplot2", "readxl", "here"))
```

These packages provide tools for plotting, data import and an easy way to
reference project files.

## Reproducing the Analysis

All of the statistical tests for the project are stored in
`scripts/middle_margin_tests.Rmd`. Knit this document to reproduce the entire
analysis:

```bash
R -e "rmarkdown::render('scripts/middle_margin_tests.Rmd')"
```

This will create an output document with the results of the tests.

## Data Directory

Raw data files are stored under `raw_data/`. The directory should contain the
following Excel spreadsheets:

* `zooplankton_count_data.xlsx` – counts of zooplankton by sample
* `zooplankton_site_data.xlsx` – metadata about the sampling sites

The repository assumes this directory structure:

```
supernova/
├── raw_data/
│   ├── zooplankton_count_data.xlsx
│   └── zooplankton_site_data.xlsx
├── scripts/
│   └── middle_margin_tests.Rmd
└── README.md
```

Ensure that the raw data files remain in `raw_data/` so that the scripts can
find them using the `here` package.
