# SSR-Spatial-Interpolation
This repository provides codes for a spatial interpolation application in a Surface Solar Radiation dataset GEBA in the paper "Machine learning techniques for spatial interpolation of solar radiation observations".

## Running the codes

### Prerequisites

You will need to install a couple of R packages before being able to run the program.

Required packages

* `tidyverse`, ` ranger`, `sp`, `raster`, `gstat`, `mgcv`, `ggplot2`, `latex2exp`, `lubridate`, `zoo`. 

Once the prerequisites above have been installed, you are ready to run the scripts.

### Program structure

The project mainly include the following:

* Four R scripts.


| File                     | Description                           |
| ------------------------ | ------------------------------------- |
| Regression_kriging_git.R | Implements regression, kriging models |
| RandomForest_git.R       | Implements random forest model        |
| modelEvaluation_git.R    | Evaluates model performance           |
| Visualization_git.R      | Plots SSR trends                      |

* The folder `data` includes a sample dataset, together with six folders corresponding to relavent model outputs/results for six continents.

* The folder `plot` is the directory to store output figures from the programs.




