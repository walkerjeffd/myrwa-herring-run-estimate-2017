Mystic River Herring Run Estimates for 2017
===========================================

Jeffrey D Walker, PhD (jeff@walkerenvres.com)
[Walker Environmental Research, LLC](https://walkerenvres.com)

March 26, 2018

## Overview

This repo contains an R script to estimate the daily and total herring run for the Mystic River
based on in-person volunteer fish counts at the Upper Mystic Lake dam.

The methodology for these estimates were obtained from MA DMF TR-25 (Nelson, 2006) using the
stratified 2-way random sampling (St2WRS) algorithm.

The results have been compared to the estimates reported by DMF using their own software, and
were found to agree exactly.

## Scripts

All calculations are performed in the `run-estimate-2017.R` script. This script requires two
packages, which can be installed by:

```r
install.packages(c("tidyverse", "lubridate"))
```

The daily and total run estimates are exported to csv files (`run-estimate-2017-daily.csv` and
`run-estimate-2017-total.csv`).
