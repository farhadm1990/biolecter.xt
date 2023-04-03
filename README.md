# Biolecter_xt R package
## An R package to generate quality plots from Biolecter XT machine.
This package is compatible with biolecter XT model output which is an excel file with different  sheets inside. The function takes directory to the excel file, a working directory for the output files, the number of sheets in the excel file (very important), and a binary (TRUE/FALSE) for the presence of biological (or technical) replicates. The funcition, then, generates timeseries plots of different filterset values over a range of speciefic time
# Installation

## 1. Isntall and library devtools pakcage on your machine
```R
install.packages("devtools")
library(devtools)
```

## 2. Download and isntall biolecter
```R
devtools::install_github("farhadm1990/biolecter.xt")
library(biolecter.xt)
```

## 3. Rinning `biolecter_xt` function
```R
biolecter.xt:::biolect_xt(working_dir=
caliberation=
n_sheets=
path_to_biolector_xlsx=
replic=
)
```

