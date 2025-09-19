assign_LCA2 <- function(x){
  # make possible LCAs
  # no family ties, assign family as LCA
  fam_LCA <- x %>% 
    pull(family) %>% 
    unique()
  genus_LCA <- x %>% 
    pull(genus) %>% 
    unique()
  species_LCA <- x %>% 
    pull(Scientific.name) %>% 
    unique()
  
  #
  if(length(fam_LCA) > 1){
    LCA <- "Uncertain"
  } else if(length(genus_LCA) > 1){
    LCA <- fam_LCA
  } else if(length(species_LCA) > 1){
    LCA <- genus_LCA
  } else {
    LCA <- species_LCA
  }
  return(LCA)
}