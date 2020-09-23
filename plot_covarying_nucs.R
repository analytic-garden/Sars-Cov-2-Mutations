plot_covarying_nucs <- function(df,  positions = c(241, 3037, 14408, 23403)) {
  require(tidyverse)
  require(rlist)
  library(grid)
  library(gridExtra)
  
  plot_list <- list()
  
  for(pos in positions) {
    col <- paste('Pos.', pos, sep='')
    df_col <- nuc_vary(df, col)
    cv_col <- df_col %>% 
      group_by_('Collection.Date', col)  %>%
      summarise(Total = sum(Count))
    
    cv_col$Collection.Date <- as.Date(cv_col$Collection.Date)
    
    title <- paste('Nucleotide', col, 'by Date')
    p <- ggplot(cv_col, aes_string(x = 'Collection.Date', y = 'Total', fill = col)) + 
      geom_bar(stat='identity', width = 1) + 
      labs(x='Collection Date', title = title)
    plot_list <- list.append(plot_list, p)
  }
  
  grid.arrange(grobs = plot_list)
}