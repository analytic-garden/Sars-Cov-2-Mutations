covary_table <- function(df) {
  require(tidyverse)
  
  df_table <- df %>% 
    select(Collection.Date, Country, Pos.241, Pos.3037, Pos.14408, Pos.23403) %>% 
    group_by(Country, Collection.Date) %>% 
    count(Pos.241, Pos.3037, Pos.14408, Pos.23403) %>% 
    mutate(Nucleotides = paste(Pos.241, Pos.3037, Pos.14408, Pos.23403,  sep = '-')) %>% 
    arrange(Collection.Date, Country) %>% 
    select(Collection.Date, Country, n, `Nucleotides`) %>% 
    rename(Count = n)
    
  return(df_table)
}