#' @title Write the pedigree file for a Wombat analysis from an R data frame.
#'
#' @description Write a pedigree file which can be used by Wombat using the
#'   information contained in an R data frame.
#'
#' @details 0 in id is considered as NA.
#'
#'   Depends on the \code{gdata} library.
#'
#' @param data a \code{data.frame} containing the data. The id should be
#'   integers. Id set to 0 are considered as NA.
#' @param animal.id,father.id,mother.id the names of the columns containing the
#'   animal, father and mother identity information (\code{strings})
#' @param file.path path to save the pedigree. Write to the console if
#'   \code{file.path} is an empty string (default empty string).
#' @param convert.pedigree \code{Boolean}, should the pedigree information be
#'   converted to multigeneration compatible pedigree using
#'   \code{\link{convert.pedigree.multigeneration}}? (default F)
#'
#' @return Nothing if \code{file.path} is given (writes to the file), writes
#'   the content of the potential pedigree file to the console if
#'   \code{file.path} is \code{""}.
#'
#' @examples
#' # generate a random pedigree structure for 20 offsprings, using 5 different
#' # sires and 7 different dams
#' set.seed(4)
#' animal_id = 1:20
#' sire_id = sample(1:5, size = 20, replace = T)
#' dam_id = sample(1:7, size = 20, replace = T)
#'
#' # generate some random traits
#' weight = rnorm(n = 20, mean = 50)
#' length = rnorm(n = 20, mean = 100)
#'
#' # assemble the information into one data frame
#' ped_pheno_data = data.frame(animal_id, sire_id, dam_id, weight, length)
#' ped_pheno_data
#'
#' # write the pedigree to the console instead of sending it to a file
#' wombat.write.pedigree(data = ped_pheno_data, animal.id = "animal_id",
#'   father.id = "sire_id", mother.id = "dam_id", convert.pedigree = T)
#'
#' @export
#'



wombat.write.pedigree = function(data, animal.id, father.id, mother.id, 
                                 file.path = "", convert.pedigree = F) {

  library(gdata) # for write.fwf at the end of the function

  # conversion if needed
  
  if (convert.pedigree) {
    
    data = convert.pedigree.multigeneration(data = data, animal.id = animal.id, 
                                            father.id = father.id, 
                                            mother.id = mother.id)
    
  }
  
  
  
  # get the pedigree information
  
  ped.data = data[, c(animal.id, father.id, mother.id)]
  
  # convert 0 to NA
  
  ped.data[ped.data == 0] = NA
  
  animals = ped.data[, animal.id]
  
  fathers = ped.data[, father.id]
  
  mothers = ped.data[, mother.id]
  
  
  
  # get the parents for which information is missing
  
  missing.fathers = unique(na.omit(fathers[!(fathers %in% animals)]))
  
  missing.mothers = unique(na.omit(mothers[!(mothers %in% animals)]))
  
  missing.parents = c(missing.fathers, missing.mothers)
  
  
  
  # create rows for those parents
  
  ped.parent = data.frame(missing.parents, 
                          rep(0, length(missing.parents)),
                          rep(0, length(missing.parents)))
  
  
  
  # concatenate with the previous pedigree
  
  ped = data.frame(ANIMAL = animals, FATHER = fathers, MOTHER = mothers)
  
  # convert NA to 0
  
  ped[is.na(ped)] = 0
  
  names(ped.parent) = names(ped)
  
  ped = rbind(ped, ped.parent)
  
  ped = ped[order(ped$ANIMAL), ]
  
  
  
  # write the pedigree (http://stackoverflow.com/questions/2851015/convert-data-frame-columns-from-factors-to-characters)
  
  ped = data.frame(lapply(ped, as.character), stringsAsFactors = F)
  
  ped = rbind(c("#ANIMAL", "FATHER", "MOTHER"), ped)
  
  write.fwf(ped, file = file.path, colnames = F)
  
}
