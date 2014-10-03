#' @title Convert a string from Wombat output to a vector
#'
#' @description Split a string from a Wombat output file into substrings based
#'   on the blank spaces.
#'
#' @details None.
#'
#' @param string a \code{string}
#'
#' @return A vector of substrings corresponding to the initial string splitted
#'   along the blank spaces.
#'
#' @examples
#' # example of Wombat output
#' output = "***** Estimates of residual covariances   ************************************\n
#'   Order of fit         =       1\n
#'   Covariance matrix\n
#'     3.4074\n
#'   Matrix of correlations and variance ratios\n
#'     0.4725\n
#'   Covariances & correlations with approximate sampling errors\n
#'     1  COVS Z 1 1        3.40737       0.549224       vrat     0.472   0.082"
#'
#' lines = strsplit(output, "\n")[[1]]
#' for (i in lines) {
#'   print(wombat.split.output.line(i))
#'  }
#'
#' @export
#'



wombat.split.output.line = function(string) {
    
  o = strsplit(string, " +")[[1]]
  
  o = o[nchar(o) > 0]
  
  o
  
}
