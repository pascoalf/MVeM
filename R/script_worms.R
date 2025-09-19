library(worms)

# fetch species data from WoRMS, dataset is "positions", name of the taxa is on the "species" column
positions$species <- gsub("sp.", "", positions$species)
positions$species <- trimws(positions$species)
unique_species <- unique(positions$species)
unique_species <- unique_species[!is.na(unique_species)]
unique_species <- unique_species[unique_species!=""]
species_worms <- wormsbymatchnames(unique_species, verbose=T, ids=F)
species_worms <- species_worms[,c("scientificname", "lsid", "rank", "kingdom", "phylum", "class", "order", "family", "genus")]
colnames(species_worms)[which(colnames(species_worms)=="scientificname")] <- "scientificName"
species_worms$species <- unique_species
positions <- join(positions, species_worms[,c("species", "scientificName")], by="species", type="left")
positions <- positions[,-which(colnames(positions)=="species")]
