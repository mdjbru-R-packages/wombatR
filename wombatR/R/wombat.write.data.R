#' @title Write the data file for a Wombat analysis from an R data frame
#'
#' @description Write a data file which can be used by Wombat using the
#'   information contained in an R data frame.
#'
#' @details 0 in id is considered as NA.
#'
#'   Depends on the \code{gdata} library.
#'
#' @param data a \code{data.frame} containing the data. The id should be
#'   integers. Id set to 0 are considered as NA.
#' @param trait the name of the column to be used to write the trait values
#'   in the data file (\code{string}).
#' @param animal.id,father.id,mother.id the names of the columns containing the
#'   animal, father and mother identity information (\code{strings})
#' @param file.path path to save the pedigree. Write to the console if
#'   \code{file.path} is an empty string (default empty string).
#' @param write.mother \code{Boolean}, should the maternal id be written also
#'   in the data file? (default T)
#' @param trait.factor multiplying factor applied to the trait before
#'   writing it (\code{float}, default 1), needed if the variance is very small
#'   otherwise Wombat will complain.
#' @param convert.pedigree \code{Boolean}, should the pedigree information be
#'   converted to multigeneration compatible pedigree using
#'   \code{\link{convert.pedigree.multigeneration}}? (default F)
#' @param filter.unknown.mother \code{Boolean}, should the individuals for
#'  which the maternal id is 0 be removed before writing the file? (default F)
#'
#' @return Nothing if \code{file.path} is given (writes to the file), writes
#'   the content of the potential data file to the console if
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
#' # send the data to the console instead of writing it to a file
#' wombat.write.data(data = ped_pheno_data, trait = "weight",
#'   animal.id = "animal_id", father.id = "sire_id", mother.id = "dam_id",
#'   write.mother = T, trait.factor = 2, convert.pedigree = T,
#'   filter.unknown.mother = T)
#'
#' @export
#' 


# Wombat write data
#------------------



wombat.write.data = function(data, trait, animal.id, father.id, mother.id, 
                             file.path = "", write.mother = T, trait.factor = 1,
                             convert.pedigree = F, filter.unknown.mother = F) {

  library(gdata)
  
  # convert the pedigree information if needed
  
  if (convert.pedigree) {
    
    data = convert.pedigree.multigeneration(data = data, animal.id = animal.id, 
                                            father.id = father.id, 
                                            mother.id = mother.id)
    
  }
  
  
  
  # prepare the data frame
  
  animals = data[, animal.id]
  
  fathers = data[, father.id]
  
  mothers = data[, mother.id]
  
  trait.values = data[, trait] * trait.factor
  
  data.dat = data.frame(ANIMAL = animals, FATHER = fathers, MOTHER = mothers,
                        VALUE_FOR_ANALYSIS = trait.values)
  
  data.dat = data.dat[!is.na(data.dat$VALUE_FOR_ANALYSIS), ] # na.omit for trait.values
    
  
  
  # format the data frame
  
  if (write.mother) {
    
    data.dat = data.dat[, c("ANIMAL", "MOTHER", "VALUE_FOR_ANALYSIS")]
    
    data.dat = data.frame(lapply(data.dat, as.character), stringsAsFactors = F)
    
    data.dat = rbind(c("#ANIMAL", "MOTHER", "VALUE_FOR_ANALYSIS"), data.dat)
    
  }
  
  else {
    
    data.dat = data.dat[, c("ANIMAL", "VALUE_FOR_ANALYSIS")]
    
    data.dat = data.frame(lapply(data.dat, as.character), stringsAsFactors = F)
    
    data.dat = rbind(c("#ANIMAL", "VALUE_FOR_ANALYSIS"), data.dat)
    
  }
  
  
  
  # write the data file
  
  if (filter.unknown.mother) {
    
    data.dat = data.dat[data.dat$MOTHER != "0", ]
    
  }
  
  write.fwf(data.dat, file = file.path, colnames = F)
  
}
