#' @title Parse the results of a permutation test analysis
#'
#' @description Parse the output of \code{\link{wombat.permutation.test.h2}}.
#'
#' @param results output from \code{\link{wombat.permutation.test.h2}}
#* @param na.rm boolean, should NAs corresponding to unsuccessful Wombat
#'   runs be removed? Default \code{FALSE}.
#'
#' @return A vector with the heritability values for each permuted dataset
#'   (including NAs if \code{na.rm=F}).
#'
#' @examples
#' # use the gryphon dataset from Wilson et al. 2010
#' #
#' head(gryphon)
#'
#' # permutations (20 perms)
#' b = wombat.permutation.test.h2(data = gryphon, trait.id = "birth.weight",
#'        animal.id = "animal", father.id = "sire", mother.id = "dam",
#'        var.animal = 1, var.residual = 1,
#'        maternal.effect = F, convert.pedigree = F,
#'        n.perm = 20)
#' c = wombat.parse.permutations.h2(b)
#' c
#'
#' @export
#'

wombat.parse.permutations.h2 = function(results, na.rm = F) {

  r = unlist(lapply(results, function(x) x[["animal"]][["vrat"]]))

  if (na.rm) {

    r = r[!is.na(r)]
    
  }
  
  return(r)
  
}
