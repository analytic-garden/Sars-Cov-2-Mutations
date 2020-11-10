plot_MI <- function(df, min_MI = 0.5, layout = 'kk') {
  # plot_MI - plot a graph of positions having significant mutual information
  # arguments:
  #   df - a data frame from MI cvs produced by MI.py
  #   min_MI - all pairs with this value or above will be used
  
  require(tidyverse)
  require(tidygraph)
  require(ggraph)
  
  temp <- df %>% 
    filter(MI >= min_MI) %>%
    mutate(from = Position_1, to = Position_2) %>%
    select(from, to, MI)
  
  mi_graph <- as_tbl_graph(temp, directed = FALSE)
  
  p <- mi_graph %>% 
        ggraph(layout = layout) + 
        geom_node_point(size = 20, color = 'steelblue') +
        geom_node_text(aes(label = name), color = 'white') +
        geom_edge_link(aes(edge_color = cut_interval(MI, 8)), edge_width = 1) +
        guides(edge_color = guide_legend(title = 'MI')) +
        ggtitle('MI Graph') +
        theme_graph(base_family = 'sans') 

  print(p)
}
  