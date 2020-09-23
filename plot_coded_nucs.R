plot_coded_nucs <- function(df, pos_list = c('Pos.3037', 'Pos.14408')) {
  require(tidyverse)
  require(ggpubr)
  require(rlist)
  
  code <- c('-' = 0, 'N' = 1, 'A' = 2, 'C' = 3, 'G' = 4, 'T' = 5)
  
  plot_list <- list()
  for(pos in pos_list) {
    df2 <- df %>% select(Collection.Date, {{ pos }})
    df2$Code <- code[pull(df2[, pos])]
    p <- ggplot(df2, aes(x = Collection.Date, y = Code)) + 
      geom_point() +
      ylim(0, 5) +
      ggtitle(pos) +
      scale_y_continuous(breaks = 0:5, 
                         labels = c('-', 'N', 'A', 'C', 'G', 'T'))
    plot_list <- list.append(plot_list, p)
  }
  
  ggarrange(plotlist = plot_list, ncol = 1)
}