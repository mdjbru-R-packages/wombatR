#' @title Parse the output file from a Wombat run (SumEstimates.out)
#'
#' @description Parse the information contained in SumEstimates.out after a
#'   Wombat run and make it available in a R object.
#'
#' @details None.
#'
#' @return A \code{list} containing the output results.
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
#' # save the files
#' wombat.write.pedigree(data = ped_pheno_data, animal.id = "animal_id",
#'   father.id = "sire_id", mother.id = "dam_id", convert.pedigree = T,
#'   file.path = "toto.ped")
#' wombat.write.data(data = ped_pheno_data, trait = "weight",
#'   animal.id = "animal_id", father.id = "sire_id", mother.id = "dam_id",
#'   write.mother = T, trait.factor = 2, convert.pedigree = T,
#'   filter.unknown.mother = T, file.path = "toto.dat")
#' wombat.write.parameters(trait = "weight", pedigree.file = "toto.ped",
#'   data.file = "toto.dat", var.residual = 2, var.animal = 2,
#'   var.mother = 2, maternal.effect = T, file.path = "toto.par")
#'
#' # run the analysis
#' wombat.unit.run(parameter.file = "toto.par")
#'
#' # parse the results
#' a = wombat.parse.result()
#* a
#'
#' @export
#' 

wombat.parse.result = function() {
  


  # check if Wombat produced a result file
  
  if (file.exists("SumEstimates.out")) {
  
    
    
    # check if there was a convergence problem
    
    f = file("SumEstimates.out", "r")
    
    content = readLines(f, -1)
    
    close(f)
    
    if (sum(unlist(lapply(strsplit(content, "Full convergence has not been achieved"), length)) > 1) == 1) {
      
      # convergence was not not achieved, return NA
      
      result = list()
      
      result$logL = NA
      
      result$residual = data.frame(variance = NA, 
                                   variance.se = NA, 
                                   vrat = NA, 
                                   vrat.se = NA)
      
      result$animal = data.frame(variance = NA, 
                                 variance.se = NA, 
                                 vrat = NA, 
                                 vrat.se = NA)
      
      result$maternal = data.frame(variance = NA, 
                                   variance.se = NA, 
                                   vrat = NA, 
                                   vrat.se = NA)
      
      result$phenotypic = data.frame(variance = NA, 
                                     variance.se = NA)
      
      result
      
    }
    
    
    
    # check if "Sum.Estimates.out" does not contain "correlations with approximate sampling"
    
    else if (!(sum(unlist(lapply(strsplit(content, "correlations with approximate sampling"), length)) > 1) >= 1)) {
      
      # results will be missing
      
      result = list()
      
      result$logL = NA
      
      result$residual = data.frame(variance = NA, 
                                   variance.se = NA, 
                                   vrat = NA, 
                                   vrat.se = NA)
      
      result$animal = data.frame(variance = NA, 
                                 variance.se = NA, 
                                 vrat = NA, 
                                 vrat.se = NA)
      
      result$maternal = data.frame(variance = NA, 
                                   variance.se = NA, 
                                   vrat = NA, 
                                   vrat.se = NA)
      
      result$phenotypic = data.frame(variance = NA, 
                                     variance.se = NA)
      
      result
            
    }
    
    
    
    else {
      
      # no convergence problem
      
      # check if a maternal effect was present
      
      f = file("SumEstimates.out", "r")
      
      content = readLines(f, -1)
      
      if (sum(unlist(lapply(strsplit(content, "Estimates for RE  2   \"MOTHER\""), length)) > 1) == 1) {
        
        maternal.effect = T
        
      }
      
      else {
        
        maternal.effect = F
        
      }
      
      close(f)
      
      # open the file
      
      f = file("SumEstimates.out", "r")
      
      # read the log likelihood
      
      found.logL = F
      
      line = readLines(f, 1)
      
      while (!found.logL) {
        
        line = readLines(f, 1)
        
        if (length(strsplit(line, "Maximum log L          =")[[1]]) > 1) {
          
          found.logL = T
          
          logL.value = line
          
        }
        
      }
      
      # read the residual value
      
      found.residual = F
      
      line = readLines(f, 1)
      
      while (!found.residual) {
        
        line = readLines(f, 1)
        
        if (length(strsplit(line, "Estimates of residual covariances")[[1]]) > 1) {
          
          found.residual = T
          
          found.value = F
          
          while (!found.value) {
            
            line = readLines(f, 1)
            
            if (length(strsplit(line, "Covariances & correlations with approximate sampling error")[[1]]) > 1) {
              
              found.value = T
              
              residual.value = readLines(f, 1)
              
            }
            
          }
          
        }
        
      }
      
      # read the ANIMAL value
      
      found.animal = F
      
      line = readLines(f, 1)
      
      while (!found.animal) {
        
        line = readLines(f, 1)
        
        if (length(strsplit(line, "Estimates for RE  1   \"ANIMAL\"")[[1]]) > 1) {
          
          found.animal = T
          
          found.value = F
          
          while (!found.value) {
            
            line = readLines(f, 1)
            
            if (length(strsplit(line, "Covariances & correlations with approximate sampling error")[[1]]) > 1) {
              
              found.value = T
              
              animal.value = readLines(f, 1)
              
            }
            
          }
          
        }
        
      }
      
      # read the maternal effect value
      
      if (maternal.effect == T) {
        
        found.maternal = F
        
        line = readLines(f, 1)
        
        while (!found.maternal) {
          
          line = readLines(f, 1)
          
          if (length(strsplit(line, "Estimates for RE  2   \"MOTHER\"")[[1]]) > 1) {
            
            found.maternal = T
            
            found.value = F
            
            while (!found.value) {
              
              line = readLines(f, 1)
              
              if (length(strsplit(line, "Covariances & correlations with approximate sampling error")[[1]]) > 1) {
                
                found.value = T
                
                maternal.value = readLines(f, 1)
                
              }
              
            }
            
          }
          
        }
        
      }
      
      # read the phenotypic value
      
      found.phenotypic = F
      
      line = readLines(f, 1)
      
      while (!found.phenotypic) {
        
        line = readLines(f, 1)
        
        if (length(strsplit(line, "Estimates of phenotypic covariances")[[1]]) > 1) {
          
          found.phenotypic = T
          
          found.value = F
          
          while (!found.value) {
            
            line = readLines(f, 1)
            
            if (length(strsplit(line, "Covariances & correlations with approximate sampling error")[[1]]) > 1) {
              
              found.value = T
              
              phenotypic.value = readLines(f, 1)
              
            }
            
          }
          
        }
        
      }
      
      # close
      
      close(f)
      
      # return
      
      result = list()
      
      result$logL = as.numeric(wombat.split.output.line(logL.value)[5])
      
      a = wombat.split.output.line(residual.value)
      
      result$residual = data.frame(variance = as.numeric(a[6]), 
                                   variance.se = as.numeric(a[7]), 
                                   vrat = as.numeric(a[9]), 
                                   vrat.se = as.numeric(a[10]))
      
      a = wombat.split.output.line(animal.value)
      
      result$animal = data.frame(variance = as.numeric(a[6]), 
                                 variance.se = as.numeric(a[7]), 
                                 vrat = as.numeric(a[9]), 
                                 vrat.se = as.numeric(a[10]))
      
      if (maternal.effect == T) {
        
        a = wombat.split.output.line(maternal.value)
        
        result$maternal = data.frame(variance = as.numeric(a[6]), 
                                     variance.se = as.numeric(a[7]), 
                                     vrat = as.numeric(a[9]), 
                                     vrat.se = as.numeric(a[10]))
        
      }
      
      a = wombat.split.output.line(phenotypic.value)
      
      result$phenotypic = data.frame(variance = as.numeric(a[6]), 
                                     variance.se = as.numeric(a[7]))
      
      result
      
    }
        
    
    
  } # end of 'if Wombat produced some results'
  
  
  
  # if Wombat did not produce any result
  
  else {
    
    result = list()
    
    result$logL = NA
    
    result$residual = data.frame(variance = NA, 
                                 variance.se = NA, 
                                 vrat = NA, 
                                 vrat.se = NA)
    
    result$animal = data.frame(variance = NA, 
                               variance.se = NA, 
                               vrat = NA, 
                               vrat.se = NA)
    
    result$maternal = data.frame(variance = NA, 
                                 variance.se = NA, 
                                 vrat = NA, 
                                 vrat.se = NA)
          
    result$phenotypic = data.frame(variance = NA, 
                                   variance.se = NA)
    
    result
    
  }
  
}

