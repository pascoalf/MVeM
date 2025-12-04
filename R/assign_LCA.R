assign_LCA <- function(x){
  # make possible LCAs
  dom_LCA <- x %>% pull(Domain) %>% unique()
  phyl_LCA <- x %>%  pull(Phylum) %>% unique()
  class_LCA <- x %>% pull(Class) %>% unique()
  ord_LCA <- x %>%  pull(Order) %>% unique()
  fam_LCA <- x %>% pull(Family) %>% unique()
  genus_LCA <- x %>%  pull(Genus) %>% unique()
  species_LCA <- x %>%pull(Species) %>% unique()
  
  taxonomic_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  #
  if(length(dom_LCA) > 1){
    LCA <- "Uncertain"
    Level <- taxonomic_levels
  } else if(length(phyl_LCA) > 1){
    LCA <- dom_LCA
    Level <- taxonomic_levels[2:7]
  } else if(length(class_LCA) > 1){
    LCA <- phyl_LCA
    Level <- taxonomic_levels[3:7]
  } else if(length(ord_LCA) > 1){
    LCA <- class_LCA
    Level <- taxonomic_levels[4:7]
  } else if(length(fam_LCA) > 1){
    LCA <- ord_LCA
    Level <- taxonomic_levels[5:7]
  } else if(length(genus_LCA) > 1){
    LCA <- fam_LCA
    Level <- taxonomic_levels[6:7]
  } else if(length(species_LCA) > 1){
    LCA <- paste(genus_LCA, "sp.")
    Level <- taxonomic_levels[7]
  } else {
    LCA <- species_LCA
    Level <- NA
  }
  return(c(LCA, Level))
}
