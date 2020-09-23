summarise.nucs <- function(cv_table) {
  require(tidyverse)
  
  temp <- cv_table %>% 
          ungroup() %>%
          filter(Nucleotides == 'C-C-C-A' | Nucleotides == 'T-T-T-G') %>%
          group_by(Collection.Date, Nucleotides) %>% 
          summarise(Total=sum(Count))
  
  return(temp)
}