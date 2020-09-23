plot_tree <- function(tree, df) {
  require(tidyverse)
  require(ggtree)
  
  p <- ggplot(tree.samp4, aes(x,y)) + geom_tree() + theme_tree() 
#  p %<+% df + geom_tiplab(aes(color=tag))
  p %<+% df + geom_tippoint(aes(color=tag, shape=usa))
}