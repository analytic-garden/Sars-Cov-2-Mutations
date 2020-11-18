covary_table2 <- function(df, positions = c(241, 3037, 14408, 23403)) {
  require(rlang)
  require(tidyverse)
  require(rlist)
  
  columns <- paste('Pos.', positions, sep = '')
  cols <- list()
  for(i in 1:length(columns)) {
    cols <- list.append(cols, sym(columns[i]))
  }

  # !! unquotes a simple argument, !!! unquotes multiple args
  df_temp <- df %>% 
    select(Collection.Date, Country, {{ columns }}) %>% 
    group_by(Country, Collection.Date) %>% 
    count(!!!cols)    # more strange NSE syntax
  
  df_temp <- df_temp %>%
    mutate(!! 'Nucleotides' := !! parse_expr(columns[1]))

  # := is needed because R doesn't allow expressions as names with '='
  for(i in 2:length(columns)) {
    df_temp <- df_temp %>% 
      mutate(!! 'Nucleotides' := paste(Nucleotides, !! parse_expr(columns[i]), sep=''))
  }
  
  df_temp <- df_temp %>%
    arrange(Collection.Date, Country) %>% 
    select(Collection.Date, Country, n, `Nucleotides`) %>% 
    rename(Count = n)
    
  return(df_temp)
}