# library(imager)
# library(shapes)
# library(StereoMorph)
# library(ggplot2)
# library(matlib)
# library(Momocs)
# library(geomorph)
# library(Morpho)
# library(abind)
# library(tibble)
# source('R/ImportSM.R')
# source('R/Functions_morpho.R')

# List of packages for session
.packages <- c(
  'imager',
  'shapes',
  'StereoMorph',
  'ggplot2',
  'matlib',
  'Momocs',
  'geomorph',
  'Morpho',
  'abind',
  'tibble',
  'intRinsic', # Intrinsic Dim
  'uwot' # UMAP
)

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session
lapply(.packages, require, character.only=TRUE)

# Load sources
source('R/ImportSM.R')
source('R/Functions_morpho.R')
