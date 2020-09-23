plot_varying_pct <- function(df) {
  # plot_varying_pct - plot Per Cent area charts for correlated genome positions
  #
  # arguments:
  # df - a dataframe produced by consensus3.py, e.g. sars_cov_2_variation_ncbi_no_dups_98.0.csv
  #
  # produce a plot of all list positions and individual .png plots for each position
  
  require(tidyverse)
  require(rlist)
  library(grid)
  library(gridExtra)
  
  plot_path = 'G:\\Covid-19\\2020_09_04\\Plots_pct\\'
  
  positions <- c('25563',
                 '1059',
                 '241',
                 '14408',
                 '23403',
                 '3037',
                 '28881',
                 '28882',
                 '28883',
                 '28144',
                 '8782',
                 '18060',
                 '17858',
                 '17747')
  
  plot_list <- list()
  
  # this code is a mess and should be simplified
  # Tidyverse doesn't do a very good job of handling columns stored in a varaible
  for(pos in positions) {
    col <- paste('Pos.', pos, sep='')
    
    col_sym = sym(col)  # WTF?
    
    tempA <- df %>% 
      select(Collection.Date, !!col_sym) %>%
      filter(!!col_sym == 'A')  %>% # WTF tidyverse syntax/semantics gets worse all the time Why .[[ ]] instead of {{ }}?
      group_by(Collection.Date) %>% 
      summarise(A = n())

    tempC <- df %>% 
      select(Collection.Date, !!col_sym) %>%
      filter(!!col_sym == 'C')  %>% 
      group_by(Collection.Date) %>% 
      summarise(C = n())
    
    tempG <- df %>% 
      select(Collection.Date, !!col_sym) %>%
      filter(!!col_sym == 'G')  %>% 
      group_by(Collection.Date) %>% 
      summarise(G = n())
    
    tempT <- df %>% 
      select(Collection.Date, !!col_sym) %>%
      filter(!!col_sym == 'T')  %>% 
      group_by(Collection.Date) %>% 
      summarise(T = n())
    
    tempN <- df %>% 
      select(Collection.Date, !!col_sym) %>%
      filter(!!col_sym == 'N')  %>% 
      group_by(Collection.Date) %>% 
      summarise(N = n())
    
    tempDel <- df %>% 
      select(Collection.Date, !!col_sym) %>%
      filter(!!col_sym == '-')  %>% 
      group_by(Collection.Date) %>% 
      summarise(Del = n())
    
    temp1 <- full_join(tempA, tempC, by = "Collection.Date") 
    temp1[is.na(temp1)] <- 0
    
    temp2 <- full_join(tempG, tempT, by = "Collection.Date") 
    temp2[is.na(temp2)] <- 0
    
    temp3 <- full_join(tempN, tempDel, by = "Collection.Date") 
    temp3[is.na(temp3)] <- 0
    
    tempj1 <- full_join(temp1, temp2, by = "Collection.Date") 
    tempj1[is.na(tempj1)] <- 0
    
    tempj2 <- full_join(tempj1, temp3, by = "Collection.Date") 
    tempj2[is.na(tempj2)] <- 0
    
    t <- tempj2 %>% pivot_longer(cols = c(A, C, G, T, N, Del), names_to = 'Nuc', values_to = "Count")
    t2 <- t %>% group_by(Collection.Date, Nuc) %>% summarise(n = sum(Count)) %>% mutate(Pct = n/sum(n))
    
    title <- paste('Nucleotide Per Cent', col)
    p <- ggplot(t2, aes(x = Collection.Date, y= Pct, fill = Nuc)) + 
      geom_area(alpha=0.6 , size=1, colour="black") + 
      labs(x='Collection Date', title = title) +
      ylim(c(0, 1))

    ggsave(paste(plot_path, 'Pos_', col, '.png', sep=''))
    plot_list <- list.append(plot_list, p)
  }
  
  grid.arrange(grobs = plot_list)
}