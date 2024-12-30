# Core tidyverse and data manipulation
library(tidyverse)  # includes dplyr, ggplot2, purrr, etc.
library(data.table)
library(dtplyr)

# File paths and project organization
library(here)

# Data import and cleaning
library(haven)      # for reading stata/sas/spss files
library(janitor)    # for data cleaning
library(ipumsr)     # for IPUMS data
library(tidycensus) # for census data

# Statistical modeling
library(msm)        # for multi-state models
library(mstate)     # for multi-state modeling
library(splines)    # for spline functions
library(survey)     # for survey analysis
library(srvyr)      # tidyverse-friendly survey analysis
library(MASS)       # for statistical functions
library(rms)        # for regression modeling
library(DEoptim)    # for optimization

# Matrix operations
library(Matrix)     # for sparse matrices
library(expm)       # for matrix exponentials

# Visualization and reporting
library(ggthemes)
library(hrbrthemes)
library(ggsci)      # scientific journal color palettes
library(directlabels)
library(patchwork)  # for combining plots
library(knitr)      # for report generation
library(kableExtra) # for table formatting
library(flextable)  # for flexible tables
library(Hmisc)      # for statistical reporting
library(glue)       # for string interpolation

# String formatting options
options(scipen = 5)

# Function aliases to avoid conflicts
transpose <- purrr::transpose
select <- dplyr::select

# Global variables
manifest = TRUE