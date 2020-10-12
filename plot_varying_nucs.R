plot_varying_nucs <- function(df, mi_table, min_MI = 0.5) {
  require(tidyverse)
  require(rlist)
  library(grid)
  library(gridExtra)
  
  plot_path = 'G:\\Covid-19\\2020_10_07\\Plots\\'
  
  temp <- mi_table %>% filter(MI >= min_MI)
  positions <- as.character(sort(unique(c(temp$Position_1, temp$Position_2))))
  
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
    ggsave(paste(plot_path, 'Pos_', col, '.png', sep=''))
    plot_list <- list.append(plot_list, p)
  }
  
  grid.arrange(grobs = plot_list)
}