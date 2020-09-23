plot_quartet <- function(cv_table, filter = FALSE) {
  require(ggplot2)
  
  title <- 'Nucleotide Quartets by Date'
  
  df <- cv_table
  if(filter) {
    df <- df %>% filter(Nucleotides == 'T-T-T-G' | Nucleotides == 'C-C-C-A')
    title <- paste(title, '- Filtered')
  }
  
  p <- ggplot(df, aes(x=as.Date(Collection.Date), y=Count, fill=Nucleotides)) + 
    geom_bar(stat='identity') + 
    labs(x='Collection Date', title = title)
         
  p
}