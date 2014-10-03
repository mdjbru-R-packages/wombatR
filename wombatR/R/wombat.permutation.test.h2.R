#' @title Permutation testing of h2 estimates
#'
#' @description Run Wombat analyses with permuted trait values.
#'
#' @details This function works by permuting the trait values among individuals
#'   and performing the Wombat analysis on permuted values. Since it happens
#'   that Wombat runs are not successful, the number of permutations can be
#'   specified in two ways, either using \code{n.perm} or
#'   \code{n.successful.perm}. If using \code{n.perm}, the number of actual
#'   heritability values can be lower than the number of permutations since
#'   it is possible that some runs do not result in a heritability estimate.
#'   When using \code{n.successful.perm}, permutations are performed until
#'   enough successful runs took place.
#'   Please note that this approach is exploratory and should be used with
#'   extreme caution: it is unclear why and when Wombat runs fail, and this
#'   might result in biased sampling of the permutation space. Hence those
#'   results should not be used to determine a significance level for
#'   heritability values as is.
#'   See also \code{\link{wombat.estimate.h2}} for additional details.
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
#' @param n.perm integer, number of permutations to be performed.
#' @param n.successful.perm integer, number of permutations resulting in
#'   successful runs to be performed.
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
#' # permutations (20 perms)
#' b = wombat.permutation.test.h2(data = gryphon, trait.id = "birth.weight",
#'        animal.id = "animal", father.id = "sire", mother.id = "dam",
#'        var.animal = 1, var.residual = 1,
#'        maternal.effect = F, convert.pedigree = F,
#'        n.perm = 20)
#' c = unlist(lapply(b, function(x) x[["animal"]][["vrat"]]))
#' c
#'
#' # idem but with 20 successful permutations
#' b = wombat.permutation.test.h2(data = gryphon, trait.id = "birth.weight",
#'        animal.id = "animal", father.id = "sire", mother.id = "dam",
#'        var.animal = 1, var.residual = 1,
#'        maternal.effect = F, convert.pedigree = F,
#'        n.successful.perm = 20)
#' c = unlist(lapply(b, function(x) x[["animal"]][["vrat"]]))
#' c
#' 
#' @export
#'

wombat.permutation.test.h2 = function(data, trait.id, 
  animal.id, father.id, mother.id,
  var.animal, var.residual, var.mother = NA,
  maternal.effect = F,
  trait.factor = 1,
  convert.pedigree = F,
  keep.tmp.directory = F,
  n.perm = NA,
  n.successful.perm = NA) {
  
  library(dplR) # for the uuid generator

  stopifnot(sum(is.na(c(n.perm, n.successful.perm))) == 1)

  if (!is.na(n.perm)) {

    output = list()

    for (i in 1:n.perm) {

      data[[trait.id]] = sample(data[[trait.id]])

      r = wombat.estimate.h2(data = data, trait.id = trait.id, 
        animal.id = animal.id, father.id = father.id, mother.id = mother.id,
        var.animal = var.animal, var.residual = var.residual,
        var.mother = var.mother, maternal.effect = maternal.effect,
        trait.factor = trait.factor,
        convert.pedigree = convert.pedigree,
        keep.tmp.directory = keep.tmp.directory)

      output[[i]] = r
      
    }

  } else {

    output = list()

    n.successful = 0

    i = 0
    
    while (n.successful < n.successful.perm) {

      i = i + 1
      
      data[[trait.id]] = sample(data[[trait.id]])

      r = wombat.estimate.h2(data = data, trait.id = trait.id, 
        animal.id = animal.id, father.id = father.id, mother.id = mother.id,
        var.animal = var.animal, var.residual = var.residual,
        var.mother = var.mother, maternal.effect = maternal.effect,
        trait.factor = trait.factor,
        convert.pedigree = convert.pedigree,
        keep.tmp.directory = keep.tmp.directory)

      output[[i]] = r
      
      value = r[["animal"]][["vrat"]]

      if (!is.na(value)) {

        n.successful = n.successful + 1
        
      }
      
    }

  }
  
  return(output)
  
}
