assign_LCA <- function(x){
  # make possible LCAs
  dom_LCA <- x %>% pull(Domain) %>% unique()
  phyl_LCA <- x %>%  pull(Phylum) %>% unique()
  class_LCA <- x %>% pull(Class) %>% unique()
  ord_LCA <- x %>%  pull(Order) %>% unique()
  fam_LCA <- x %>% pull(Family) %>% unique()
  genus_LCA <- x %>%  pull(Genus) %>% unique()
  species_LCA <- x %>%pull(Species) %>% unique()
  
  #
  if(length(dom_LCA) > 1){
    LCA <- "Uncertain"
    Level <- "Domain"
  } else if(length(phyl_LCA) > 1){
    LCA <- dom_LCA
    Level <- "Phylum"
  } else if(length(class_LCA) > 1){
    LCA <- phyl_LCA
    Level <- "Class"
  } else if(length(ord_LCA) > 1){
    LCA <- class_LCA
    Level <- "Order"
  } else if(length(fam_LCA) > 1){
    LCA <- ord_LCA
    Level <- "Family"
  } else if(length(genus_LCA) > 1){
    LCA <- fam_LCA
    Level <- "Genus"
  } else if(length(species_LCA) > 1){
    LCA <- paste(genus_LCA, "sp.")
    Level <- "Species"
  } else {
    LCA <- species_LCA
    Level <- NA
  }
  return(c(LCA, Level))
}
