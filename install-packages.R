packages <- c("targets","tarchetypes", "conflicted",
              "terra", "tidyverse", "dismo", "extRemes",
              "digest", "fitdistrplus", "devtools",
              "future", "future.callr", "clustermq",
              "parallel", "RColorBrewer", "gbm",
              "tidyterra", "ggpubr", "PMCMRplus",
              "xtable", "ggVennDiagram", "ggcorrplot",
              "splancs", "ggh4x", "ggokabeito", "ggnewscale",
              "wesanderson")

lapply(packages, function(somepackage){
if(!require(somepackage, character.only = T)){install.packages(somepackage)}
})

devtools::install_github("sjevelazco/flexsdm@HEAD")
