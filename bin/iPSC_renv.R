#### This script can help you to setup RENV for this project to manage packages across (virtual) machines. ####

#### 1) Make sure you updated R to the current working version (R4.4.0 as of 14-05-2024)
#### 2) Add the renv library (found on the Z:\software\R) to your R library of the current working version (C:\Program Files\R\R-4.4.0\library)
#### 3) Activate renv on your (V)M for this version of R, see 'setup renv' below
#### 4) Check if your project is up-to-date and keep it updated using the provided functions in 'keep your project updated'



# setup renv
library(renv)
renv::init() # make a renv project (already done for this project, so no longer needed)
renv::activate() # Only use once on a (V)M to setup renv for that version of R

# keep your project updated
renv::status() # Check if your packages are up-to-date compared to the project renv.lock file
renv::restore() # restore missing packages (only works when you can install packages from within the DRE again, otherwise add the packages manually to the renv library folder in your project)
renv::snapshot() # Use this to add packages that you have installed to the renv.lock file of the project
renv::install() # install a new packages that you require, instead of the install.package() function

# other helpful functions
renv::dependencies() # Will check which packages are required for all the scripts within the whole project
renv::update() # update all packages (only works when you can install packages from within the DRE again, otherwise add the packages manually to the renv library folder in your project)
renv::remove() # opposite of renv::install(), to remove a package that is not used within the whole project