# Check ties function
check_ties <- function(x){
  x %>%
    pull(Species) %>% 
    unique() %>% 
    word() %>%
    unique() 
} 
