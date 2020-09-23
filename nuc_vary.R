nuc_vary <- function(df, column) {
  # nuc_vary - produce a count of each nucleotide by country and date
  # 
  # arguments:
  # df - a dataframe produced by consensus3.py, e.g. sars_cov_2_variation_ncbi_no_dups_98.0.csv
  # column - a genome column number in format Pos.xxxx
  #
  # returns:
  # df_table - a dataframe of Country Collection.Date, Count, column 
  #            Count is a count of the total counts of a nucleotide for that date and Country
  #
  # tidyverse is a mess when it comes to hadling clumn names as variables
  # NSE syntax is inconsistent. Why does count() operate differently?
  
  require(tidyverse)
  
  col_sym = sym(column)
  df_table <- df %>% 
    select('Collection.Date', 'Country', !!col_sym) %>% 
    group_by(Country, Collection.Date) %>% 
    count(!!col_sym) %>%     # WTF Have to use deprecated function because count({{ column }}) add quotes to existing column name!
    arrange(Collection.Date, Country) %>% 
    select('Collection.Date', 'Country', 'n', !!col_sym) %>% 
    rename(Count = n)
  
  return(df_table)
}