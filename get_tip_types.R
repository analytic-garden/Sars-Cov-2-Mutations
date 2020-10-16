get_tip_types <- function(tree) {
  require(tidyverse)
  
  tips <- tree$tip.label
  df <- data.frame(seq = tips, 
                   tag = as.factor(str_extract(tips, '.-.-.-.')),
                   usa = ! is.na(str_extract(tips, 'USA')))
  
  return(df)
}