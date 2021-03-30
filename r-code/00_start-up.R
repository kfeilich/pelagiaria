# 00_start-up.R

# Author: Kara Feilich
# Date modified: 26 November 2018
# This script loads necessary packages and functions from external sources
# for subsequent analysis.

# NB: spptest is available here: https://rdrr.io/github/myllym/spptest/f/README.md
needed.packages <- list('ape', 'phytools', 'geiger', 'here', 'purrr', 'spptest',
                        'dplyr', 'ggplot2', 'plyr', 'grid', 'lattice', 
                        'gridExtra', 'ggplotify', 'stringr', 'rfishbase', 
                        'readr','tidyr')


# NEED to use CORRECTED source code from https://github.com/mwpennell/geiger-v2/blob/master/R/disparity.R
source("https://raw.githubusercontent.com/mwpennell/geiger-v2/master/R/disparity.R")
# Modified dtt function to call modified MDI function and allow user to change y-axis
source("https://raw.githubusercontent.com/djmurrell/DTT-Envelope-code/master/dtt1.R")
# Modified code for two sided MDI test. Note this is the CORRECTED version of the MDI test in geiger
source("https://raw.githubusercontent.com/djmurrell/DTT-Envelope-code/master/getMDI1.R")
# This function takes a dtt produced object and runs it through the rank envelope test
source("https://raw.githubusercontent.com/djmurrell/DTT-Envelope-code/master/rank_dtt.R")

# For plotting timescale
source("https://gist.githubusercontent.com/willgearty/fd46c71ee8335f47ac95a4552f34715c/raw/7704de3fc8947cc3cf529edfa362a6b2d2032d40/gggeo_scale.R")

lapply(needed.packages, library, character.only = TRUE)
rm(needed.packages)  # clean up

theme_custom <- function() {
  theme_bw(base_size = 12, base_family = "TimesNewRoman")  %+replace% 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
}



