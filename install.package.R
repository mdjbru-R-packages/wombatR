# creating a package
# http://adv-r.had.co.nz/Package-quick-reference.html
# note: the latest version of devtools should be installed after installing
# the CRAN one by using
# library(devtools)
# install_github("devtools")
# to avoid a bug about "keep.source"
#
# to be able to run the script with Rscript, see
# http://stackoverflow.com/questions/8964515/cannot-call-roxygenize-function-from-rscript-batch-file


#library(devtools)

#create("mtdnar") # this creates a directory for the package files

#create("matti.toolkit") # this creates a directory for the package files


#dev_mode() # to install development packages aside from the regular ones

#load_all("mtdnar") # to load the development package



# function to update the package under development

package.names = commandArgs(trailingOnly = TRUE)

update.dev.package = function(package.name) {
  
  library(devtools)

  library(methods)

  library(utils)
  
  # package.name should also be the path to the package, i.e. the folder
  # containing the package should be in the working directory (but not
  # the working directory itself)
  
  try(library(package.name, character.only = T))
  
  try(unload(package.name))
  
  try(remove.packages(package.name))
  
  document(package.name)
  
  install(package.name)
  
  library(package.name, character.only = T)
  
}



# install

a = sapply(package.names, update.dev.package)



