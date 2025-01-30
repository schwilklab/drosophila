# Data and code for "Spermless males fail to remedy the low fecundity of parthenogenetic females in Drosophila melanogaster"

Authors: Lewis I. Held, Jr.*, Surya J. Banerjee, Dylan W. Schwilk, Souvik Roy, Kambre A. Huddleston, and Jason J. Shin.

Manuscript resubmitted February 2025.


## To run analysis
  - data are in csv files in the /data directory
  - two R script files are in the /scripts directory
  - run `stats.R` to produce all statistical output (to console) and the manuscript figure (saved as both pdf and png file in the results directory).
  
## Dependencies

The code relies on the following R packages: 

```R
library(readr)
library(tidyr)
library(dplyr)
library(survival)
library(coxme) # for mixed effect Cox model
library(glmmTMB)  # for poisson regression.
library(survminer) # for prettier plots using ggplot2 - not necessary
library(car) 
library(patchwork) # for figure
```
