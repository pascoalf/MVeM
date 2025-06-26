# function to check ties
check_ties <- function(x){
  x %>%
    pull(Scientific.name) %>% 
    unique() %>% 
    word() %>%
    unique() 
}