# R-FDA Project
This project focuses on implementing and modifying various statistical methods from the [R Views Statistics](https://rviews.rstudio.com/categories/statistics/) tutorials using CDISC datasets.

## Overview
Welcome to the **R-FDA Project**! This repository contains R scripts for functional data analysis, survival analysis, longitudinal data analysis, and advanced visualizations using CDISC datasets. The project aims to adapt methods from the R Views Statistics tutorials to clinical data programming.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Directory Structure](#directory-structure)
- [Implemented Methods](#implemented-methods)
- [Contributing](#contributing)
- [License](#license)

## Prerequisites
Ensure you have the following installed:
- R (version 4.0 or higher)
- RStudio (recommended)
- Required R packages:
  - `tidyverse`, `devtools`, `FDA`, `survival`, `lme4`, `haven`, `ggplot2`, `dplyr`, `rmarkdown`, `shiny`, `DT`, `lubridate`

You can install the packages using:

```r
install.packages(c('tidyverse', 'devtools', 'FDA', 'survival', 'lme4', 'haven', 'ggplot2', 'dplyr', 'rmarkdown', 'shiny', 'DT', 'lubridate'))
```

## Installation
Clone the repository to your local machine:

```bash
git clone https://github.com/Sakshi16-code/R-FDA_Project.git
cd R-FDA_Project
```

Set your working directory in RStudio:

```r
setwd('path/to/R-FDA_Project')
```

## Usage
### Example Scripts from R Views:
- Functional Data Analysis: `scripts/FDA/fda_analysis.R`
- Survival Analysis: `scripts/Survival_Analysis/cox_model.R`
- Longitudinal Analysis: `scripts/Longitudinal/mixed_model.R`
- Advanced Visualizations: `scripts/Visualization/ggplot_custom.R`

Run a script in R by sourcing it:

```r
source('scripts/Survival_Analysis/cox_model.R')
```

## Directory Structure
```bash
R-FDA_Project/
├── outputs/                  # Generated output files (reports, plots, tables)
├── scripts/                 
│   ├── FDA/                 # Functional Data Analysis
│   │   └── fda_analysis.R
│   ├── Survival_Analysis/   # Survival models (Kaplan-Meier, Cox PH)
│   │   └── cox_model.R
│   ├── Longitudinal/        # Mixed-effects models, time series
│   │   └── mixed_model.R
│   └── Visualization/       # ggplot2 visualizations
│       └── ggplot_custom.R
├── README.md                # Project documentation
└── requirements.txt         # Required R packages
```

## Implemented Methods
### Functional Data Analysis
- Smoothing and basis expansion
- Functional regression models

### Survival Analysis
- Kaplan-Meier estimates
- Cox proportional hazards models

### Longitudinal Data Analysis
- Mixed-effects models
- Time series decomposition

### Visualization
- Advanced visualizations using ggplot2
- R Markdown reports for reproducible analysis

## Contributing
Contributions are welcome! Follow these steps to contribute:
1. Fork the repository.
2. Create a new branch: `git checkout -b feature-branch`
3. Make your changes and commit: `git commit -m 'Add new feature'`
4. Push to your branch: `git push origin feature-branch`
5. Open a pull request.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements
- [R Views](https://rviews.rstudio.com/) for statistical methods and tutorials.
- [CDISC](https://www.cdisc.org/) for clinical data standards.
- [RStudio](https://www.rstudio.com/) for supporting the R community.
