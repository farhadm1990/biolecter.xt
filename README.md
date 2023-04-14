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
biolecter.xt:::biolecter_xt(path_to_biolector_xlsx="path_to_file.xlsx",
working_dir="output_directory",
caliberation= "both", #to show calibrated (TRUE), raw (FALSE) or both of them "both".
n_sheets= 8, #Number of sheets in the output xlsx dataset (oblicatory).
replic= TRUE #if you have biological/technical replica
)
```

### An example output plot

![fig1](https://github.com/farhadm1990/biolecter.xt/blob/master/output/p_Timeseries_plot_of_calibrated%20and%20raw_data_for_isolate_Js15-Glucose-pH_with_5_filtersets%20(Anaerobic%20mode%20On).%20Replicate_date:Wed%20Apr%2012%2023:14:30%202023.jpeg)

