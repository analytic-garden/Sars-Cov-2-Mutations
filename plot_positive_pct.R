plot_positive_pct <- function(lag = 7, counties = c('Albany', 'Columbia', 'Rensselaer', 'Saratoga', 'Schenectady')) {
  require(tidyverse)
  require(zoo)
  require(rlist)
  require(grid)
  require(gridExtra)
  
  ny_url <- "https://health.data.ny.gov/api/views/xdss-u53e/rows.csv?accessType=DOWNLOAD&api_foundry=true"
  df <- read.csv(ny_url)
  df$Test.Date <- as.Date(df$Test.Date, format='%m/%d/%Y')
  
  plot_list <- list()
  df2 <- data.frame(Test.Date = NULL, County = NULL, Pct = NULL, Avg = NULL)
  for(in_county in counties) {
    temp <- df %>% 
      filter(County == {{ in_county }}) %>% 
      mutate(Pct = ifelse(Total.Number.of.Tests.Performed == 0, 0, (New.Positives/Total.Number.of.Tests.Performed) * 100)) %>%
      select(Test.Date, County, Pct) %>%
      mutate(Avg = rollapply(Pct, lag, mean, align='right', fill = NA)) %>%
      mutate(Avg = ifelse(Avg >= 0, Avg, 0))
    
    p2 <- ggplot(temp, aes(x=Test.Date, y=Avg)) + 
          geom_line() +
          labs(title = paste('Pct Positive Cases', {{ in_county }}, lag, 'Day Moving Average', sep = ' '),
               y = 'Avg')
    plot_list <- list.append(plot_list, p2)
    
    df2 <- rbind(df2, temp)
  }

  p <- ggplot(df2, aes(x=Test.Date, y=Avg, color = County)) + 
      geom_line() +
      labs(title = paste('Pct Positive Cases', lag, 'Day Moving Average', sep = ' '),
           y = 'Avg')
  print(p)
  
  grid.arrange(grobs = plot_list, ncol=2, top = textGrob(Sys.Date()))
  
  return(df2)
}