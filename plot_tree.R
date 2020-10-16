plot_tree <- function(tree, df) {
  require(tidyverse)
  require(ggtree)
  
  # %<+% is a ggtree operator to add annotation data to a  tree
  
  p <- ggplot(tree, aes(x,y)) + geom_tree() + theme_tree() 
#  p %<+% df + geom_tiplab(aes(color=tag))
  p %<+% df + geom_tippoint(aes(color=tag, shape=usa))
}