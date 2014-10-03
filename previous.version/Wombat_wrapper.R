# source for PKD trout Estonia 2011
# (works with R version 2.15.3 (2013-03-01) and
# Wombat Version 19-03-2013)



library(gdata) # to write tables in fixed width format (write.fwf)

library(dplR) # to have an uuid generator



#-----------------------------------------------------------------------------#
#                        NON WOMBAT-SPECIFIC FUNCTIONS                        #
#-----------------------------------------------------------------------------#



# convert pedigree with independent parent and offspring id numbering to 
# multigeneration compatible one
#-----------------------------------------------------------------------



convert.pedigree.multigeneration = function(data, animal.id, father.id, 
                                            mother.id) {
  
  # Convert the id for animal, father and mother from an independent 
  # format (i.e. id start at 1 for animal, father and mother) to
  # a multigeneration compatible format (i.e. no overlap between
  # offspring and parent id).
  #
  # Note: 0 in id is considered as NA
  #
  # IMPORTANT: This should NOT be used on multigeneration pedigrees 
  # where offspring and parents id are related, since such relations
  # will be broken during the recoding of the id.
  #
  # TAKES
  # data: a data frame
  # animal.id, father.id, mother.id: the names of the columns containing the
  #   animal, father and mother information
  #
  # RETURNS
  # the same data frame with updated id
  
  
  
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



#-----------------------------------------------------------------------------#
#                     LOW-LEVEL WOMBAT-RELATED FUNCTIONS                      #
#-----------------------------------------------------------------------------#



# convert a string from Wombat output to a vector
#------------------------------------------------



wombat.split.output.line = function(string) {
  
  # Split a string from Wombat output file
  #
  # TAKES
  # string : output string from Wombat
  #
  # RETURNS
  # a vector of strings with the elements from string splitted by blank spaces
  
  
  
  o = strsplit(string, " +")[[1]]
  
  o = o[nchar(o) > 0]
  
  o
  
}



# Wombat write pedigree
#----------------------



wombat.write.pedigree = function(data, animal.id, father.id, mother.id, 
                                 file.path = "", convert.pedigree = F) {
  
  # Write the pedigree file for Wombat analysis from a data frame.
  #
  # A 0 in the id is considered as NA.
  #
  # Depends on library(gdata).
  #
  # TAKES
  # data: a data frame
  # animal.id, father.id, mother.id: the names of the columns containing the
  #   animal, father and mother information
  # file.path: a path to save the pedigree, write to the console if empty string.
  # convert.pedigree: boolean. Should the pedigree information be converted to
  #   multigeneration compatible pedigree using 
  #   convert.pedigree.multigeneration() ?
  #
  # RETURNS
  # nothing if file.path is given (writes to the file), writes to the console
  # the content of the potential pedigree file is file.path is "".
  
  
  
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



# Wombat write data
#------------------



wombat.write.data = function(data, trait, animal.id, father.id, mother.id, 
                             file.path = "", write.mother = T, trait.factor = 1,
                             convert.pedigree = F, filter.unknown.mother = F) {
  
  # Write the data file for Wombat analysis from a data frame.
  # A 0 in the id is considered as NA.
  # 
  # Depends on library(gdata).
  #
  # TAKES
  # data: a data frame
  # trait: name of the column containing the trait of interest
  # animal.id, father.id, mother.id: the names of the columns containing the
  #   animal, father and mother information
  # file.path: a path to save the data, write to the console if empty string.
  # write.mother: boolean, should it write the mother id to the data file?
  # trait.factor: numerical, multiplying factor applied to the trait column, 
  #   needed if the variance is very small otherwise Wombat will complain.
  # convert.pedigree: boolean. Should the pedigree information be converted to
  #   multigeneration compatible pedigree using 
  #   convert.pedigree.multigeneration() ?
  # filter.unknown.mother: boolean, if TRUE removes the rows with MOTHER 0.
  #
  # RETURNS
  # nothing if file.path is given (writes to the file), writes to the console
  # the content of the potential data file is file.path is "".
  


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



# generate a parameter file for Wombat (univariate) #
#---------------------------------------------------#



wombat.write.parameters = function(trait, pedigree.file, data.file, var.residual,
                                 var.animal, var.mother = NA, var.random.effect = NA,
                                 maternal.effect = F,
                                 random.effect = F, file.path) {
  
  # Produce a parameter input file for Wombat
  #
  # TAKES
  # trait : string, trait to be analysed
  # pedigree.file : path to the pedigree file
  # data.file : path to the data file
  # var.residual : initial value for residual variance
  # var.animal : initial value for animal variance
  # var.mother : initial value for mother effect variance
  # maternal.effect : boolean to include mother as a random factor
  # random.effect : boolean to include a random.effect (must also have been included
  #   with wombat.data.save)
  # file.path : path to save the parameter file
  
  
  
  # open the file
  
  f = file(file.path, "w")
  
  cat("# A single optional comment line of up to 74 characters\n", file = f)
  
  cat("COMMENT WOMBAT ", trait, "\n", file = f, sep = "")
  
  cat("# The type of analysis to be performed. In this case a UNIvariate analysis\n", file = f)
  
  cat("ANALYSIS UNI\n", file = f)
  
  cat("# Name the pedigree file, assuming it is in the same folder as the parameter file\n", file = f)
  
  cat("PEDS ", pedigree.file, "\n", file = f, sep = "")
  
  cat("# Name the data file, assuming it is in the same folder as the parameter file\n", file = f)
  
  cat("DATA ", data.file, "\n", file = f, sep = "")
  
  cat("# Description of the data file\n", file = f)
  
  cat("# Each variable to be fitted in the model needs to be followed by\n", file = f)
  
  cat("# the maximum number of levels\n", file = f)
  
  cat("ANIMAL 1000\n", file = f)
  
  cat("MOTHER 1000\n", file = f)
  
  if (random.effect) {
    
    cat("RANDOM_EFFECT_TRAIT\n", file = f)
    
  }
  
  cat("VALUE_FOR_ANALYSIS", "\n", file = f, sep = "")
  
  cat("END DATA\n", file = f)
  
  cat("# Model specification\n", file = f)
  
  cat("# Type of effect (FIXed, COVariate, RANdom) and variable name\n", file = f)
  
  cat("# NRM indicates that a pedigree is available for ANIMAL\n", file = f)
  
  cat("MODEL\n", file = f)
  
  cat("RAN ANIMAL NRM\n", file = f)
  
  if (random.effect == T) {
    
    cat("RAN RANDOM_EFFECT_TRAIT\n", file = f)
    
  }
  
  if (maternal.effect == T) {
    
    cat("RAN MOTHER\n", file = f)
    
  }
  
  cat("# The trait to be analysed\n", file = f)
  
  cat("TRAIT ", "VALUE_FOR_ANALYSIS", "\n", file = f, sep = "")
  
  cat("END\n", file = f)
  
  cat("# Specify starting values and the number of rows and columns in the matrix\n", file = f)
  
  cat("# These are not very important here\n", file = f)
  
  cat("VAR residual 1\n", file = f)
  
  cat(var.residual, "\n", file = f, sep = "")
  
  cat("VAR ANIMAL 1\n", file = f)
  
  cat(var.animal, "\n", file = f, sep = "")
  
  if (random.effect == T) {
    
    cat("VAR RANDOM_EFFECT_TRAIT 1\n", file = f)
    
    cat(var.random.effect, "\n", file = f, sep = "")
    
  }
  
  if (maternal.effect == T) {
    
    cat("VAR MOTHER 1\n", file = f)
    
    cat(var.mother, "\n", file = f, sep = "")
    
  }
  
  close(f)
  
}



# run Wombat for one elementary analysis
#---------------------------------------



wombat.unit.run = function(parameter.file) {
  
  # Run Wombat
  #
  # TAKES
  # parameter.file : the path to the parameter file
  
  
  
  command_line = paste("wombat", parameter.file, sep = " ")
  
  system(command_line, intern = T)
  
}



# read wombat output file
#------------------------



wombat.parse.result = function() {
  
  # Parse the output file from Wombat (SumEstimates.out)
  #
  # TAKES
  # nothing
  #
  # RETURNS
  # a list containing the information
  
  
  
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



#-----------------------------------------------------------------------------#
#                     HIGH-LEVEL WOMBAT-RELATED FUNCTIONS                     #
#-----------------------------------------------------------------------------#



wombat.estimate.heritability = function(data, trait.id, 
                                        animal.id, father.id, mother.id,
                                        var.animal, var.residual, var.mother = NA,
                                        maternal.effect = F,
                                        trait.factor = 1,
                                        convert.pedigree = F) {
  
  # Perform a Wombat run to estimate heritability for one trait from a data 
  # frame containing the trait values and the pedigree information
  #
  # TAKES
  # data: a data frame
  # trait.id: column name for the trait to be analysed
  # animal.id, father.id, mother.id: the names of the columns containing the
  #   animal, father and mother information
  # var.animal: initial value for animal variance
  # var.residual: initial value for residual variance
  # var.mother: initial value for mother effect variance
  # maternal.effect: boolean to include mother as a random factor
  # trait.factor: numerical, multiplying factor applied to the trait column, 
  #   needed if the variance is very small otherwise Wombat will complain.
  # convert.pedigree: boolean. Should the pedigree information be converted to
  #   multigeneration compatible pedigree using 
  #   convert.pedigree.multigeneration() ?
  
  
  
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
                          var.random.effect = NA, maternal.effect = maternal.effect, 
                          random.effect = F, file.path = f.par)
  
  
  
  # run Wombat
  
  setwd(working.dir)
  
  wombat.unit.run(parameter.file = basename(f.par))
  
  
  
  # parse the results
  
  results = wombat.parse.result()
  
  
  
  # delete the working directory
  
  setwd("..")
  
  unlink(working.dir, recursive = T)
  
  
  
  # return
  
  results
  
}



# Manual tests
#-------------



if (F) {
  
  d = read.table("../data/tables/master.table", header = T, sep = "\t")
  
  trait = "CAT"
  
  wombat.write.pedigree(d, "fish_global_id", "col1_parent_0", "col1_parent_1", 
                        "toto.ped", T)
  
  wombat.write.data(d, trait,  "fish_global_id", "col1_parent_0", "col1_parent_1", 
                    "toto.dat", trait.factor = 1, convert.pedigree = T)
  
  wombat.write.parameters(trait, "toto.ped", "toto.dat", var.residual = 250,
                          var.animal = 38, file.path = "toto.par")
  
  wombat.unit.run("toto.par")

  #a = wombat.parse.result()
  
  #a = wombat.estimate.heritability(d, "swollenness", "fish_global_id", 
  #                                 "col1_parent_0","col1_parent_1",30, 200, 
  #                                 NA, F, 100, T)
  
  # effect of prior on results
  
  trait = "swollenness"
  
  trait.factor = 100
  
  trait.var = var(d[, trait] * trait.factor, na.rm = T)

  var.animal.frac = seq(0.1, 0.9, by = 0.1)
 
  var.animal = var.animal.frac * trait.var

  var.residuals = trait.var - var.animal

  h2 = vector()
  
  h2.se = vector()
  
  for (i in 1:length(var.animal.frac)) {
    
    a = wombat.estimate.heritability(d, trait, "fish_global_id", 
                                     "col1_parent_0","col1_parent_1", 
                                     var.animal[i], var.residuals[i], 
                                     var.residuals[i], T, trait.factor, T)
    
    h2 = c(h2, a$animal[3])
    
    h2.se = c(h2.se, a$animal[4])
    
  }
  
  plot(var.animal.frac, h2, xlim = c(0, 1), ylim = c(0, 1), bty = "n", pch = 16)
  
  abline(lm(as.vector(h2, mode = "numeric") ~ as.vector(var.animal.frac)))
  
}



# Tests using the gryphon dataset
#--------------------------------



# Data files and tutorial are downloaded from http://www.wildanimalmodels.org/tiki-index.php?page=the%20ecologists%20guide%20to%20the%20animal%20model#WOMBAT
# NOTE : the space character between "#" and "ANIMAL" in the gryphon_uni.dat 
# has to be removed !
# The aim is to reproduce the results from the tutorial using the R wrapper.



if (F) {

  # load the data
  
  ped = read.table("gryphon.ped_WOMBAT", header = T, comment.char = "")
  
  dat = read.table("gryphon_uni.dat", header = T, comment.char = "")
  
  
  
  # merge the data to a single data frame with trait and pedigree information
  
  data = merge(ped, dat, by.x = "X.ANIMAL", by.y = "X.ANIMAL", all.x = T, all.y = T)
  
  
  
  # estimate heritability of birth weight
  
  a = wombat.estimate.heritability(data = data, 
                                   trait.id = "BWT", 
                                   animal.id = "X.ANIMAL", 
                                   father.id = "FATHER", 
                                   mother.id = "MOTHER.x", 
                                   var.animal = 1, 
                                   var.residual = 1,
                                   var.mother = NA, 
                                   maternal.effect = F, 
                                   trait.factor = 1,
                                   convert.pedigree = F)
  
  # with a maternal effect
  
  # estimate heritability of birth weight
  
  a = wombat.estimate.heritability(data = data, 
                                   trait.id = "BWT", 
                                   animal.id = "X.ANIMAL", 
                                   father.id = "FATHER", 
                                   mother.id = "MOTHER.x", 
                                   var.animal = 1, 
                                   var.residual = 1,
                                   var.mother = 1, 
                                   maternal.effect = T, 
                                   trait.factor = 1,
                                   convert.pedigree = F)
  
}
