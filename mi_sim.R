mi_sim <- function(df, mi_table, columns = c(28881, 28882), num_samples = 1000, target = NULL) {
  # mi_sim - permutation test of significance of MI between two columns
  #
  # arguments:
  #   df - a data frame produced by variation.py
  #   mi_table - a data frame of mutual information between columns produced by MI.py
  #   columns - a pair of columns 
  #   num_samples - number of permutation tests to perform
  #   target - if NULL, use MI from mi_table other wise a number
  #
  # returns:
  #   a list
  #   prob - estimate of probability of obtaining an MI >= target from permutation test
  #   mi - a vector of MI calculations from permutation tests
  #
  # requires:
  #   target > 0 if target is not NULL
  #
  require(infotheo)
  require(tidyverse)
  
  N = dim(df)[1]
  
  if(is.null(target)) {
    temp = mi_table %>% filter(Position_1 == columns[1] & Position_2 == columns[2])
    if(dim(temp)[1] == 0) {
      temp = mi_table %>% filter(Position_1 == columns[2] & Position_2 == columns[1])
    }
    target = temp$MI
  }
  
  col1 = df[, paste('Pos.', columns[1], sep = '')]
  col2 = df[, paste('Pos.', columns[2], sep = '')]
  
  mi = rep(0, num_samples)
  c = 0
  for(i in 1:num_samples) {
    c1 = sample(col1, N, replace = TRUE)
    c2 = sample(col2, N, replace = TRUE)
    mi[i] = natstobits(mutinformation(c1, c2))
    if(mi[i] >= target) {
      c = c + 1
    }
  }
  
  return(list(prob = c / N, mi_samples = mi))
}