plot_quartet_nucs_pct <- function(df) {
  require(tidyverse)
  require(rlist)
  library(grid)
  library(gridExtra)
  
  quartet = c('Pos.241', 'Pos.3037', 'Pos.14408', 'Pos.23403')
  plot_list <- list()
  
  for(col in quartet) {
     df_col <- nuc_vary(df, col)
     cv_col <- df_col %>% 
       group_by_('Collection.Date', col)  %>%
       summarise(Total = sum(Count))
     
     cv_col$Collection.Date <- as.Date(cv_col$Collection.Date)
     
     title <- paste('Nucleotide', col, 'by Date')
     p <- ggplot(cv_col, aes_string(x = 'Collection.Date', y = 'Total', fill = col)) + 
       geom_bar(stat='identity') + 
       labs(x='Collection Date', title = title)
     plot_list <- list.append(plot_list, p)
  }
  
  grid.arrange(grobs = plot_list)
}