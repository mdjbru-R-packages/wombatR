#' @title Test the sensitivity of h2 estimates to starting values
#'
#' @description Run Wombat analyses with several starting variance values in
#'   order to test the sensitivity of h2 estimates to starting conditions.
#'
#' @details See \code{\link{wombat.estimate.h2}}.
#'
#' @param data a \code{data.frame} containing the data.
#' @param trait.id the name of the column containing the trait (\code{string})
#' @param animal.id,father.id,mother.id strings, the names of the columns
#'   containing the animal, father and mother identities (id coded as integers).
#' @param var.animal numerical vector, initial values for animal variance.
#' @param var.residual numerical vector, initial values for residual variance.
#' @param var.mother numerical vector, initial values for maternal effect
#'   variance.
#' @param maternal.effect \code{boolean}, should maternal effect be included
#'   as a random factor?
#' @param trait.factor multiplying factor applied to the trait before
#'   writing it (\code{float}, default 1), needed if the variance is very small
#'   otherwise Wombat will complain.
#' @param convert.pedigree \code{boolean}, should the pedigree information be
#'   converted to multigeneration compatible pedigree using
#'   \code{\link{convert.pedigree.multigeneration}}? (default F)
#' @param keep.tmp.directory \code{boolean}, should the temporary directories
#'   used for the analysis be kept instead of being deleted? If True, the
#'   return value will also contain the name of this directory.
#'
#' @return A list of lists containing the Wombat run results (i.e. a list
#'   of outputs from \code{link{wombat.estimate.h2}}.
#'
#' @examples
#' # use the gryphon dataset from Wilson et al. 2010
#' #
#' head(gryphon)
#' a = wombat.estimate.h2(data = gryphon, trait.id = "birth.weight",
#'        animal.id = "animal", father.id = "sire", mother.id = "dam",
#'        var.animal = 1, var.residual = 1, maternal.effect = F,
#'        convert.pedigree = F)
#' # the numerical values should mirror those obtained in Wilson's tutorial.
#'
#' # sensitivity test
#' var.obs = var(gryphon[["birth.weight"]], na.rm = T)
#' var.obs
#' var.animal = seq(1, 7, by = 1)
#' var.residual = var.obs - var.animal
#' b = wombat.sensitivity.h2(data = gryphon, trait.id = "birth.weight",
#'        animal.id = "animal", father.id = "sire", mother.id = "dam",
#'        var.animal = var.animal, var.residual = var.residual,
#'        maternal.effect = F, convert.pedigree = F)
#' c = lapply(b, function(x) c(var.animal.start = x[["var.animal.start"]],
#'                             h2 = x[["animal"]][["vrat"]]))
#' c
#'
#' @export
#'

wombat.sensitivity.h2 = function(data, trait.id, 
  animal.id, father.id, mother.id,
  var.animal, var.residual, var.mother = NA,
  maternal.effect = F,
  trait.factor = 1,
  convert.pedigree = F,
  keep.tmp.directory = F) {
  
  library(dplR) # for the uuid generator

  stopifnot(length(var.animal) == length(var.residual))

  if (!is.na(var.mother)) {
    
    stopifnot(length(var.animal) == length(var.mother))

  } else {

    var.mother = rep(NA, length(var.animal))
    
  }
  
  output = list()

  for (i in 1:length(var.animal)) {

    va = var.animal[i]

    vr = var.residual[i]

    vm = var.mother[i]

    r = wombat.estimate.h2(data = data, trait.id = trait.id, 
      animal.id = animal.id, father.id = father.id, mother.id = mother.id,
      var.animal = va, var.residual = vr,
      var.mother = vm, maternal.effect = maternal.effect,
      trait.factor = trait.factor,
      convert.pedigree = convert.pedigree,
      keep.tmp.directory = keep.tmp.directory)

    output[[i]] = r
          
  }

  return(output)
  
}
