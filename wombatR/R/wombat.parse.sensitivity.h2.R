#' @title Parse the results of a sensitivity analysis
#'
#' @description Parse the output of \code{\link{wombat.sensitivity.h2}}.
#'
#' @param results output from \code{\link{wombat.sensitivity.h2}}
#'
#' @return A data frame with two columns: the starting variance values and 
#'   estimated h2 values
#'
#' @examples
#' # use the gryphon dataset from Wilson et al. 2010
#' #
#' head(gryphon)

#' # sensitivity test
#' var.obs = var(gryphon[["birth.weight"]], na.rm = T)
#' var.obs
#' var.animal = seq(1, 7, by = 1)
#' var.residual = var.obs - var.animal
#' b = wombat.sensitivity.h2(data = gryphon, trait.id = "birth.weight",
#'        animal.id = "animal", father.id = "sire", mother.id = "dam",
#'        var.animal = var.animal, var.residual = var.residual,
#'        maternal.effect = F, convert.pedigree = F)
#' c = wombat.parse.sensitivity.h2(b)
#' c
#'
#' @export
#'

wombat.parse.sensitivity.h2 = function(results) {

  var.animal = unlist(lapply(results, function(x) x[["var.animal.start"]]))
  h2 = unlist(lapply(results, function(x) x[["animal"]][["vrat"]]))
  return(data.frame(var.animal, h2))
  
}
