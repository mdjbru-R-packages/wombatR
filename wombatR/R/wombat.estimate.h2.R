#' @title Estimate heritability with Wombat
#'
#' @description Use Wombat (REML approach) to estimate a trait heritability.
#'   Maternal effect can be included in the analysis.
#'
#' @details This is a simple R interface to call Wombat from R scripts. Wombat
#'   has to be installed on the system (see
#'   http://didgeridoo.une.edu.au/km/wombat.php).
#'
#'   For traits with very small numerical values, a multiplying
#'   \code{trait.factor} has to be applied for numerical stability.
#'
#'   When giving initial values for residuals, those values should not
#'   be greater than the observed variance of the trait (after
#'   multiplication by \code{trait.factor}).
#'
#'   The run is performed in a temporary folder whose name is generated
#'   by an uuid generator from the \code{dplR} library. It is deleted
#'   at the end of the analysis.
#'
#' @param data a \code{data.frame} containing the data.
#' @param trait.id the name of the column containing the trait (\code{string})
#' @param animal.id,father.id,mother.id strings, the names of the columns
#'   containing the animal, father and mother identities (id coded as integers).
#' @param var.animal numerical, initial value for animal variance.
#' @param var.residual numerical, initial value for residual variance.
#' @param var.mother numerical, initial value for maternal effect variance.
#' @param maternal.effect \code{boolean}, should maternal effect be included
#'   as a random factor?
#' @param trait.factor multiplying factor applied to the trait before
#'   writing it (\code{float}, default 1), needed if the variance is very small
#'   otherwise Wombat will complain.
#' @param convert.pedigree \code{boolean}, should the pedigree information be
#'   converted to multigeneration compatible pedigree using
#'   \code{\link{convert.pedigree.multigeneration}}? (default F)
#' @param keep.tmp.directory \code{boolean}, should the temporary directory
#'   used for the analysis be kept instead of being deleted? If True, the
#'   return value will also contain the name of this directory.
#'
#' @return A list containing the Wombat run results.
#'
#' @examples
#' # generate a random pedigree structure for 20 offsprings, using 5 different
#' # sires and 7 different dams
#' #
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
#' # run Wombat
#' a = wombat.estimate.h2(data = ped_pheno_data, trait.id = "weight",
#'       animal.id = "animal_id", father.id = "sire_id",
#'       mother.id = "dam_id", var.animal = 1, var.residual = 1,
#'       var.mother = 1, maternal.effect = T, trait.factor = 2,
#'       convert.pedigree = T)
#'
#' # use the gryphon dataset from Wilson et al. 2010
#' #
#' head(gryphon)
#' b = wombat.estimate.h2(data = gryphon, trait.id = "birth.weight",
#'        animal.id = "animal", father.id = "sire", mother.id = "dam",
#'        var.animal = 1, var.residual = 1, maternal.effect = F,
#'        convert.pedigree = F)
#' # the numerical values should mirror those obtained in Wilson's tutorial.
#'
#' # with maternal effect
#' c = wombat.estimate.h2(data = gryphon, trait.id = "birth.weight",
#'        animal.id = "animal", father.id = "sire", mother.id = "dam",
#'        var.animal = 1, var.residual = 1, var.mother = 1,
#'        maternal.effect = T,
#'        convert.pedigree = F)
#'
#' @export
#'

wombat.estimate.h2 = function(data, trait.id, 
  animal.id, father.id, mother.id,
  var.animal, var.residual, var.mother = NA,
  maternal.effect = F,
  trait.factor = 1,
  convert.pedigree = F,
  keep.tmp.directory = F) {
  
  library(dplR) # for the uuid generator
  
  # create a working directory
  
  directory.suffix = uuid.gen()()
  
  working.dir = paste("wombat.", directory.suffix, sep = "")
  
  dir.create(working.dir, showWarnings = F)
  
  
  
  # write the pedigree file
  
  f.ped = file.path(working.dir, "run.ped")
  
  wombat.write.pedigree(data = data, animal.id = animal.id, 
                        father.id = father.id, mother.id = mother.id, 
                        file.path = f.ped,
                        convert.pedigree = convert.pedigree)
  
  
  
  # write the data file
  
  f.dat = file.path(working.dir, "run.dat")
  
  if (maternal.effect) {
    
    filter.unknown.mother = T
    
  }
  
  else {
    
    filter.unknown.mother = F
    
  }
  
  wombat.write.data(data = data, trait = trait.id, animal.id = animal.id, 
                    father.id = father.id, mother.id = mother.id, 
                    file.path = f.dat, write.mother = T, 
                    trait.factor = trait.factor,
                    convert.pedigree = convert.pedigree,
                    filter.unknown.mother = filter.unknown.mother)
  
  
  
  # write the parameter file
  
  f.par = file.path(working.dir, "run.par")
  
  wombat.write.parameters(trait = trait.id, pedigree.file = basename(f.ped), 
                          data.file = basename(f.dat), var.residual = var.residual, 
                          var.animal = var.animal, var.mother = var.mother, 
                          maternal.effect = maternal.effect, file.path = f.par)
  
  
  
  # run Wombat
  
  setwd(working.dir)
  
  wombat.unit.run(parameter.file = basename(f.par))
  
  
  
  # parse the results
  
  results = wombat.parse.result()
  
  
  
  # delete the working directory
  
  setwd("..")

  if (!keep.tmp.directory) {
    
    unlink(working.dir, recursive = T)
  
  } else {

    results[["tmp.dir"]] = working.dir
    
  }


  
  # return

  results[["var.animal.start"]] = var.animal

  results[["var.residual.start"]] = var.residual

  results[["var.mother.start"]] = var.mother
  
  return(results)
  
}
