plot_varying_nucs2 <- function(df, positions = c(241, 3037, 14408, 23403)) {
  # plot_varying_nucs2 - plot bar charts for correlated genome positions
  #
  # arguments:
  # df - a dataframe produced by consensus3.py, e.g. sars_cov_2_variation_ncbi_no_dups_98.0.csv
  #
  # produce a plot of all list positions and individual .png plots for each position
  
  require(tidyverse)
  require(rlist)
  library(grid)
  library(gridExtra)
  library(rlist)
  library(rlang)
  
  plot_list <- list()
  
  # this code is a mess and should be simplified
  # Tidyverse doesn't do a very good job of handling columns stored in a variable
  for(pos in positions) {
    col <- paste('Pos.', pos, sep='')
    
    col_sym = sym(col)  # WTF?
    temp_list = list()
    for(let in c('A', 'C', 'G', 'T', 'N')) {
      temp_list[[let]] <- df %>% 
        select(Collection.Date, all_of(col)) %>%
#        select(Collection.Date, !!col_sym) %>%
#        filter(!!col_sym == let)  %>% # WTF tidyverse syntax/semantics gets worse all the time.
        filter(get(col) == let)  %>% # WTF tidyverse syntax/semantics gets worse all the time.
        group_by(Collection.Date) %>% 
        summarise(!!let := n())   # more strange syntax
    }
    
    # handle Indels separately
    temp_list[['Indel']] <- df %>% 
      select(Collection.Date, !!col_sym) %>%
      filter(!!col_sym == '-')  %>% 
      group_by(Collection.Date) %>% 
      summarise(Indel = n())
    
    temp1 <- full_join(temp_list[['A']], temp_list[['C']], by = "Collection.Date") 
    temp1[is.na(temp1)] <- 0
    
    temp2 <- full_join(temp_list[['G']], temp_list[['T']], by = "Collection.Date") 
    temp2[is.na(temp2)] <- 0
    
    temp3 <- full_join(temp_list[['N']], temp_list[['Indel']], by = "Collection.Date") 
    temp3[is.na(temp3)] <- 0
    
    tempj1 <- full_join(temp1, temp2, by = "Collection.Date") 
    tempj1[is.na(tempj1)] <- 0
    
    tempj2 <- full_join(tempj1, temp3, by = "Collection.Date") 
    tempj2[is.na(tempj2)] <- 0
    
    t <- tempj2 %>% pivot_longer(cols = c(A, C, G, T, N, Indel), names_to = "Nuc", values_to = "Count")
    t2 <- t %>% group_by(Collection.Date, Nuc) %>% summarise(Total = sum(Count))
    
    title <- paste('Nucleotide Counts', col)
    p <- ggplot(t2, aes(x = Collection.Date, y = Total, fill = Nuc)) + 
      geom_bar(stat='identity', width = 1) + 
      labs(x='Collection Date', title = title)

    plot_list <- list.append(plot_list, p)
  }
  
  grid.arrange(grobs = plot_list)
}