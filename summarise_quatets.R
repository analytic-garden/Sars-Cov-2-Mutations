summarise.quartets <- function(cv_table) {
  require(tidyverse)
  
  temp <- cv_table %>% 
          ungroup() %>%
          group_by(Nucleotides) %>% 
          summarise(Total=sum(Count))
  
  return(temp)
}