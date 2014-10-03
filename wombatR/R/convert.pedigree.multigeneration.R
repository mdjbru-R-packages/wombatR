#' @title Convert parents and offspring identities to unique identities
#'
#' @description Convert a pedigree with independent numbering for the identities
#'   of parents and offsprings to a multigeneration compatible
#'   pedigree, where identical numbers are not shared between different
#'   individuals between the offspring and parent columns.
#'
#' @details Convert the id for animal, father and mother from an independent 
#'   format (i.e. id start at 1 for animal, father and mother) to
#'   a multigeneration compatible format (i.e. no overlap between
#'   offspring and parent identities).
#'
#'   Note: 0 in id is considered as NA.
#'
#'   IMPORTANT: This should NOT be used on multigeneration pedigrees 
#'   where offspring and parents id are related, since such relations
#'   will be broken during the recoding of the id.
#'
#' @param data a \code{data.frame} containing the data. The id should be
#'   integers. Id set to 0 are considered as NA.
#' @param animal.id,father.id,mother.id the names of the columns containing the
#'   animal, father and mother identity information (\code{strings})
#'
#' @return A \code{data.frame} with the same data with updated id.
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
#' # convert the pedigree part of the data frame
#' ped_pheno_data2 = convert.pedigree.multigeneration(data = ped_pheno_data,
#'   animal.id = "animal_id", father.id = "sire_id", mother.id = "dam_id")
#' ped_pheno_data2
#'
#' @export
#' 



convert.pedigree.multigeneration = function(data, animal.id, father.id, 
                                            mother.id) {
  
  # get the ids
  
  data.id = data[, c(animal.id, father.id, mother.id)]
  
  # convert 0 to NA in the id
  
  data.id[data.id == 0] = NA
  
  animals = data.id[, animal.id]
  
  fathers = data.id[, father.id]
  
  mothers = data.id[, mother.id]
  
  
  
  # determine the number of individuals within each category
  
  n.animals = max(animals, na.rm = T)
  
  n.fathers = max(fathers, na.rm = T)
  
  n.mothers = max(mothers, na.rm = T)
  
  
  
  # recode the ids
  
  fathers = fathers
  
  mothers = mothers + n.fathers
  
  animals = animals + n.fathers + n.mothers
  
  ped = data.frame(animals = animals, fathers = fathers, mothers = mothers)
  
  ped[is.na(ped)] = 0
  
  
  
  # replace the columns in the original data frame
  
  data[, animal.id] = ped$animals
  
  data[, father.id] = ped$fathers
  
  data[, mother.id] = ped$mothers
  
  
  
  # return
  
  data
  
}
