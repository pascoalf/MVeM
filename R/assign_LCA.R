# function to assign LCA  
assign_LCA <- function(x){
  # identify ASV with ties to break
  ties <- check_ties(x)
  
  if(length(ties) == 1){
    # untie within the same genus
    LCA = paste(ties, "sp.")
    
  } else {
    # check if ties are from Delphinidae family
    if(mean(ties %in% delphinidae_family$Genus) == 1){
      LCA = "Delphinidae sp."
    } else if(mean(ties %in% pleuronectidae_family$Genus) == 1){
      LCA = "Pleuronectidae sp."
    } else if(mean(ties %in% ziphiidae_family$Genus)){
      LCA = "Ziphiidae sp."
    } else if(mean(ties %in% salmonidae_family$Genus)){
    LCA = "Salmonidae sp."
    } else if(mean(ties %in% mugilidae_family$Genus)){
      LCA = "Mugilidae sp."
    } else {
      LCA = "Uncertain"
    }
  }
  return(LCA)
}