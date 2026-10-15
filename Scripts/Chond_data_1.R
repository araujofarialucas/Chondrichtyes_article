# Load Packages
library(dplyr)
library(stringr)
library(rfishbase)
library(readxl)
library(ape)
library(phytools)
library(sf)
library(randtip)
library(picante)

# Load Files

## Open IUCN data
dat_iucn <- read.csv("Data/Primary data/iucn_data.csv") %>%
  dplyr::select(scientificName, redlistCategory) %>%
  rename(Species = scientificName) %>%
  rename(IUCN = redlistCategory)

## Open IUCN data taxonomy
dat_iucn_taxonomy <- read.csv("Data/Primary data/iucn_taxonomy.csv", sep = ",") %>%
  dplyr::select(scientificName, orderName, familyName, genusName) %>%
  rename(Species = scientificName) %>%
  rename(Order = orderName) %>%
  rename(Family = familyName)

## Load instrinsic traits
fisb_data <- species(
  dat_iucn$Species, 
  fields = c("Species", "Length", "DepthRangeShallow", "DepthRangeDeep")
)

fisb_data <- fisb_data %>%
  na.omit(Length, DepthRangeShallow, DepthRangeDeep) %>% 
  dplyr::select(
    Species = Species,
    Length = Length,
    `Upper Depth` = DepthRangeShallow,
    `Lower Depth` = DepthRangeDeep
  )

iucn_eco <- dat_iucn[dat_iucn$Species %in% 
                       setdiff(dat_iucn$Species, fisb_data$Species), ] %>%
  dplyr::select(Species)

write.csv(iucn_eco, "Data/Secondary data/Manual inspection/manual_inspection_eco.csv")

iucn_eco <- read_excel("Data/Secondary data/Manual inspection/iucn_eco.xlsx") %>%
  dplyr::select(Species, Length, `Upper Depth`, `Lower Depth`) %>%
  na.omit()
iucn_eco$`Lower Depth` <- as.numeric(iucn_eco$`Lower Depth`)
iucn_eco$`Upper Depth` <- as.numeric(iucn_eco$`Upper Depth`)
iucn_eco$Length <- as.numeric(iucn_eco$Length)

eco_data <- bind_rows(iucn_eco, fisb_data)

## Load Phylogeny
shark_trees <- read.nexus("Data/Trees/10.cal.tree.nex")
set.seed(123)
trees_IDs <- sample(length(shark_trees), 100)
shark_trees_final <- shark_trees[trees_IDs]
shark_trees_final <- lapply(shark_trees_final, force.ultrametric, method = "nnls")

### Creating a vector for phylogenies' species
phylo_species <- as.data.frame(shark_trees_final[[1]]$tip.label) %>%
  rename(Species = `shark_trees_final[[1]]$tip.label`)
phylo_species$Species <- gsub("_", " ", phylo_species$Species)

## Load shapefiles
distrib <- st_read("Data/Primary data/Chond_distrib") %>%
  rename(Species = SCI_NAME)

distrib <- distrib[match(unique(distrib$Species), distrib$Species), ]

## Open IUCN DD data
dat_iucn_dd <- read.csv("Data/Primary data/DD_data/dat_iucn_dd.csv") %>%
  dplyr::select(scientificName, redlistCategory) %>%
  rename(Species = scientificName) %>%
  rename(IUCN = redlistCategory)

## Open IUCN DD data taxonomy
dat_iucn_taxonomy_dd <- read.csv("Data/Primary data/DD_data/taxonomy_dd.csv", sep = ",") %>%
  dplyr::select(scientificName, orderName, familyName, genusName, authority) %>%
  rename(Species = scientificName) %>%
  rename(Ordem = orderName) %>%
  rename(Familia = familyName) %>%
  rename(Genero = genusName) %>%
  rename(Author = authority)

## Load instrinsic traits for DD
fisb_data_dd <- species(
  dat_iucn_dd$Species, 
  fields = c("Species", "Length", "DepthRangeShallow", "DepthRangeDeep")
)

fisb_data_dd <- fisb_data_dd %>%
  na.omit(Length, DepthRangeShallow, DepthRangeDeep) %>% 
  dplyr::select(
    Species = Species,
    Length = Length,
    `Upper Depth` = DepthRangeShallow,
    `Lower Depth` = DepthRangeDeep
  )

iucn_eco_dd <- dat_iucn_dd[dat_iucn_dd$Species %in% 
                             setdiff(dat_iucn_dd$Species, fisb_data_dd$Species), ] %>%
  dplyr::select(Species)

write.csv(iucn_eco_dd, "Data/Secondary data/Manual inspection/test_dd.csv")

iucn_eco_dd <- read_excel("Data/Secondary data/Manual inspection/iucn_eco_dd.xlsx") %>%
  dplyr::select(Species, Length, `Upper Depth`, `Lower Depth`) %>%
  na.omit()
iucn_eco_dd$`Lower Depth` <- as.numeric(iucn_eco_dd$`Lower Depth`)
iucn_eco_dd$`Upper Depth` <- as.numeric(iucn_eco_dd$`Upper Depth`)
iucn_eco_dd$Length <- as.numeric(iucn_eco_dd$Length)

eco_data_dd <- bind_rows(iucn_eco_dd, fisb_data_dd)

## Load shapefiles for DD
distrib_dd <- st_read("Data/Primary data/DD_distribution") %>%
  rename(Species = SCI_NAME)

distrib_dd <- distrib_dd[match(unique(distrib_dd$Species), distrib_dd$Species), ]

##################################################################################################
# Aligning IUCN data taxonomy

## Loading valid names
dat_iucn$fishbase = tryCatch(validate_names(dat_iucn$Species), error = function(e) {
  sapply(dat_iucn$Species, function(sp)
  {tryCatch(rfishbase::validate_names(sp), error = function(e) NA)})})

## Looking at the valid names with more than one correspondent in the IUCN database and the NA values for valid names
dat_iucn_error <- dat_iucn %>% 
  add_count(fishbase, name = "count") %>% 
  filter(count > 1 | is.na(fishbase)) %>%
  dplyr::select(Species, IUCN, fishbase)

write.csv(dat_iucn_error, "/Applications/Chondrichtyes_article/Data/Taxonomy_problems/iucn_problem.csv")

### Removing invalid species (same valid name than other valid species)
remove_species_iucn <- c("Harriotta haeckeli", "Squalus acutipinnis")
dat_iucn <- dat_iucn[!dat_iucn$Species %in% remove_species_iucn, ]

### Removing species not found in the fishbase by manual inspection
dat_iucn <- dat_iucn[dat_iucn$Species != "Ctenacis magnificum", ]

## Looking at species with names to alter
dat_iucn_change <- dat_iucn %>%
  dplyr::select(Species, fishbase) %>%
  filter(Species != fishbase) ##19 to alter

## Dataset with valid names
dat_iucn <- dat_iucn %>%
  mutate(Species = if_else(Species != fishbase, fishbase, Species),
         Species = gsub(" ", "_", Species))

##################################################################################################
# Aligning iucn taxonomy data taxonomy

## Loading valid names
dat_iucn_taxonomy$fishbase <- tryCatch(rfishbase::validate_names(dat_iucn_taxonomy$Species), error = function(e) {
  sapply(dat_iucn_taxonomy$Species, function(sp)
  {tryCatch(rfishbase::validate_names(sp), error = function(e) NA)})})

### Removing invalid species (same valid name than other valid species)
dat_iucn_taxonomy <- dat_iucn_taxonomy[!dat_iucn_taxonomy$Species %in% remove_species_iucn, ]

### Removing species not found in the fishbase by manual inspection
dat_iucn_taxonomy <- dat_iucn_taxonomy[dat_iucn_taxonomy$Species !=
                                         "Ctenacis magnificum", ]  #### Não foi achado

## Dataset with valid names
dat_iucn_taxonomy <- dat_iucn_taxonomy %>%
  mutate(Species = if_else(Species != fishbase, fishbase, Species),
         Species = gsub(" ", "_", Species))

##################################################################################################
# Aligning intrinsic data taxonomy

## Loading valid names
eco_data$fishbase <- tryCatch(rfishbase::validate_names(eco_data$Species), error = function(e) {
  sapply(eco_data$Species, function(sp)
  {tryCatch(rfishbase::validate_names(sp), error = function(e) NA)})})

### Removing invalid species (same valid name than other valid species)
eco_data <- eco_data[!eco_data$Species %in% remove_species_iucn, ]

### Removing species not found in the fishbase by manual inspection
eco_data <- eco_data[eco_data$Species != "Ctenacis magnificum", ]  #### Não foi achado

## Dataset with valid names
eco_data <- eco_data %>%
  mutate(Species = if_else(Species != fishbase, fishbase, Species),
         Species = gsub(" ", "_", Species))

##################################################################################################
# Aligning shp taxonomy

## Loading valid names
distrib$fishbase = tryCatch(validate_names(distrib$Species), error = function(e) {
  sapply(distrib$Species, function(sp)
  {tryCatch(rfishbase::validate_names(sp), error = function(e) NA)})})

## Looking at the valid names with more than one correspondent in the IUCN database and the NA values for valid names
distrib_prob <- distrib %>% 
  st_set_geometry(NULL) %>%
  add_count(fishbase, name = "count") %>% 
  filter(count > 1 | is.na(fishbase)) %>%
  dplyr::select(Species, fishbase)

write.csv(distrib_prob, "/Applications/Chondrichtyes_article/Data/Taxonomy_problems/distrib_problem.csv")

### Removing invalid species (same valid name than other valid species)
remove_species_distrib <- c("Harriotta haeckeli", "Squalus acutipinnis")
distrib <- distrib[!distrib$Species %in% remove_species_distrib, ]

### Removing species not found in the fishbase by manual inspection
distrib <- distrib[distrib$Species != "Ctenacis magnificum", ]  #### Não foi achado

## Looking at species with names to alter
distrib_change <- distrib %>%
  dplyr::select(Species, fishbase) %>%
  filter(Species != fishbase) ##19 to alter

## Dataset with valid names
distrib <- distrib %>%
  mutate(Species = if_else(Species != fishbase, fishbase, Species),
         Species = gsub(" ", "_", Species))

##################################################################################################
# Aligning phylogeny species taxonomy

## Loading valid names
phylo_species$fishbase = tryCatch(validate_names(phylo_species$Species), error = function(e) {
  sapply(phylo_species$Species, function(sp)
  {tryCatch(rfishbase::validate_names(sp), error = function(e) NA)})})

## Looking at the valid names with more than one correspondent in the IUCN database and the NA values for valid names
dat_phylo_prob <- phylo_species %>% 
  add_count(fishbase, name = "count") %>% 
  filter(count > 1 | is.na(fishbase)) %>%
  select(Species, fishbase)

write.csv(dat_phylo_prob, "/Applications/Chondrichtyes_article/Data/Taxonomy_problems/filog_problem.csv")

### Removing invalid species (same valid name than other valid species)
remove_species_phylo <- c("Harriotta haeckeli", "Hydrolagus lemures", "Leucoraja virginica",
                          "Leucoraja caribbaea", "Amblyraja robertsi", "Amblyraja badia",
                          "Raja africana", "Raja atriventralis", "Narcine firma", "Torpedo zugmayeri",
                          "Torpedo polleni", "Trygonorrhina melaleuca", "Rhinoptera adspersa",
                          "Mobula rochebrunei", "Mobula japanica", "Mobula diabolus", "Potamotrygon dumerilii",
                          "Potamotrygon castexi", "Dasyatis matsubarai",
                          "Dasyatis ushiei", "Dasyatis thetidis", "Proscymnodon plunketi",
                          "Isistius labialis", "Cephaloscyllium pardelotum", "Cephaloscyllium maculatum",
                          "Scyliorhinus tokubee", "Parmaturus sp.", "Mustelus walkeri"
)
phylo_species <- phylo_species [!phylo_species $Species %in% remove_species_phylo, ]

### Replacing false NAs (valid names found by manual inspection at fishbase website)
position <- which(phylo_species$Species == "Raja macrorynchus")
phylo_species[position, 2] <- "Dipturus batis"

position <- which(phylo_species$Species == "Torpedo alexandrinsis")
phylo_species[position, 2] <- "Torpedo alexandrinsis"

position <- which(phylo_species$Species == "Rhinoptera peli")
phylo_species[position, 2] <- "Rhinoptera marginata"

position <- which(phylo_species$Species == "Rhinoptera sewelli")
phylo_species[position, 2] <- "Rhinoptera jayakari"

position <- which(phylo_species$Species == "Proscyllium venustum")
phylo_species[position, 2] <- "Proscyllium venustum"

position <- which(phylo_species$Species == "Carcharhinus wheeleri")
phylo_species[position, 2] <- "Carcharhinus amblyrhynchos"

### Removing species not found in the fishbase by manual inspection
true_na_tree <- c("Rajella alia", "Okamejei porosa", "Zearaja argentinensis",
                  "Dipturus intermedia", "Narcine nigra", "Narcine bicolor",
                  "Myliobatis rhombus", "Rhinoptera hainanensis",
                  "Etmopterus schmidti", "Glyphis sp.1")

phylo_species <- phylo_species[!phylo_species$Species %in% true_na_tree, ]

### Removing invalid species (same valid name than other valid species) after changing NAs for real valid names
dat_phylo_prob <- phylo_species %>% 
  add_count(fishbase, name = "count") %>% 
  filter(count > 1 | is.na(fishbase)) %>%
  select(Species, fishbase)

remove_species_phylo <- c("Raja macrorynchus", "Rhinoptera peli", "Rhinoptera sewelli",
                          "Carcharhinus wheeleri")
phylo_species <- phylo_species [!phylo_species $Species %in% remove_species_phylo, ]

## Looking at names to alter
phylo_species_change <- phylo_species %>%
  select(Species, fishbase) %>%
  filter(Species != fishbase) ##148 to alter

## Phylogeny with correct taxonomy

### Altering incorrect names on phylogenies species' list
phylo_species_clean <- phylo_species %>%
  mutate(
    Species_valid = if_else(Species != fishbase, fishbase, Species),
    Species = gsub(" ", "_", Species),
    Species_valid = gsub(" ", "_", Species_valid)
  )

### Removing invalid species from the set of trees and changing incorrect anmes
shark_trees_clean <- lapply(shark_trees_final, function(tree) {
  remover <- setdiff(tree$tip.label, phylo_species_clean$Species)
  tree <- drop.tip(tree, remover)
  index <- match(tree$tip.label, phylo_species_clean$Species)
  tree$tip.label <- phylo_species_clean$Species_valid[index]
  return(tree)
})

###Converting the object from a list to a multiPhylo
shark_trees_clean <- do.call(c, shark_trees_clean)

##################################################################################################
# Aligning IUCN DD data taxonomy

## Loading valid names
dat_iucn_dd$fishbase = tryCatch(validate_names(dat_iucn_dd$Species), error = function(e) {
  sapply(dat_iucn_dd$Species, function(sp)
  {tryCatch(rfishbase::validate_names(sp), error = function(e) NA)})})

## Looking at the valid names with more than one correspondent in the IUCN database and the NA values for valid names
dat_iucn_dd_prob <- dat_iucn_dd %>% 
  add_count(fishbase, name = "count") %>% 
  filter(count > 1 | is.na(fishbase)) %>%
  dplyr::select(Species, IUCN, fishbase)

## Looking at species with names to alter
dat_iucn_dd_change <- dat_iucn_dd %>%
  dplyr::select(Species, fishbase) %>%
  filter(Species != fishbase) ##1 to alter

## Dataset with valid names
dat_iucn_dd <- dat_iucn_dd %>%
  mutate(Species = if_else(Species != fishbase, fishbase, Species),
         Species = gsub(" ", "_", Species))

##################################################################################################
# Aligning iucn DD taxonomy data taxonomy

## Loading valid names
dat_iucn_taxonomy_dd$fishbase <- tryCatch(rfishbase::validate_names(dat_iucn_taxonomy_dd$Species), error = function(e) {
  sapply(dat_iucn_taxonomy_dd$Species, function(sp)
  {tryCatch(rfishbase::validate_names(sp), error = function(e) NA)})})

## Dataset with valid names
dat_iucn_taxonomy_dd <- dat_iucn_taxonomy_dd %>%
  mutate(Species = if_else(Species != fishbase, fishbase, Species),
         Species = gsub(" ", "_", Species))

##################################################################################################
# Aligning intrinsic data taxonomy for DD

## Loading valid names
eco_data_dd$fishbase <- tryCatch(rfishbase::validate_names(eco_data_dd$Species), error = function(e) {
  sapply(eco_data_dd$Species, function(sp)
  {tryCatch(rfishbase::validate_names(sp), error = function(e) NA)})})

## Dataset with valid names
eco_data_dd <- eco_data_dd %>%
  mutate(Species = if_else(Species != fishbase, fishbase, Species),
         Species = gsub(" ", "_", Species))

##################################################################################################
# Aligning shp taxonomy for DD

## Loading valid names
distrib_dd$fishbase = tryCatch(validate_names(distrib_dd$Species), error = function(e) {
  sapply(distrib_dd$Species, function(sp)
  {tryCatch(rfishbase::validate_names(sp), error = function(e) NA)})})

## Looking at the valid names with more than one correspondent in the IUCN database and the NA values for valid names
distrib_dd_prob <- distrib_dd %>% 
  st_set_geometry(NULL) %>%
  add_count(fishbase, name = "count") %>% 
  filter(count > 1 | is.na(fishbase)) %>%
  dplyr::select(Species, fishbase)

## Looking at species with names to alter
distrib_dd_change <- distrib_dd %>%
  dplyr::select(Species, fishbase) %>%
  filter(Species != fishbase) ##1 to alter

## Dataset with valid names
distrib_dd <- distrib_dd %>%
  mutate(Species = if_else(Species != fishbase, fishbase, Species),
         Species = gsub(" ", "_", Species))

##################################################################################################

## Load CAS Species
taxa_completa <- load_taxa()

chondrichthyes_lista <- taxa_completa %>%
  filter(SuperClass == "Chondrichthyes") %>%
  dplyr::select(Genus, Species, SpecCode, Family)

iucn_total <- gsub("_", " ", c(dat_iucn_dd$Species, dat_iucn$Species))

new_chond <- setdiff(chondrichthyes_lista$Species, iucn_total) %>%
  gsub(pattern = "_", replacement = " ", x = .)

## Load CAS intrinsic data

fisb_data_cas <- rfishbase::species(
  new_chond, 
  fields = c("Species", "Length", "DepthRangeShallow", "DepthRangeDeep")
)

fisb_data_cas <- fisb_data_cas %>%
  na.omit(Length, DepthRangeShallow, DepthRangeDeep) %>% ##108 espécies
  dplyr::select(
    Species = Species,
    Length = Length,
    `Upper Depth` = DepthRangeShallow,
    `Lower Depth` = DepthRangeDeep
  ) %>%
  mutate(Species = gsub(" ", "_", Species))

##################################################################################################
# New shp

iucn_common <- Reduce(intersect, list(dat_iucn$Species, dat_iucn_taxonomy$Species,
                                      eco_data$Species))

iucn_dd_common <- Reduce(intersect, list(dat_iucn_dd$Species, dat_iucn_taxonomy_dd$Species,
                                         eco_data_dd$Species))

total_dd <- c(iucn_dd_common, fisb_data_cas$Species) ##131 DD species

total <- c(iucn_common, iucn_dd_common, fisb_data_cas$Species)

distrib_total <- c(distrib_dd$Species, distrib$Species)

new_shp <- setdiff(total, distrib_total)

write.csv(new_shp, "Data/Secondary data/Manual inspection/Species without shp/species_without_shp.csv")

new_shp <- st_read("Data/Primary data/New_shp") %>%
  rename(Species = binomial)

new_shp <- new_shp[match(unique(new_shp$Species), new_shp$Species), ] %>%
  mutate(Species = gsub(" ", "_", Species))

## New shps for each group
new_shp_iucn <- intersect(new_shp$Species, dat_iucn$Species)
new_shp_iucn_dd <- intersect(new_shp$Species, dat_iucn_dd$Species)
new_shp_cas <- intersect(new_shp$Species, fisb_data_cas$Species)

##################################################################################################
# Phylogeny imputation

## IUCN Species with intrinsic traits, IUCN dara, and distribution available
distrib_sp_iucn <- c(new_shp_iucn, distrib$Species)
distrib_sp_iucn_dd <- c(new_shp_iucn_dd, distrib_dd$Species)

iucn_total <- Reduce(intersect, list(dat_iucn$Species, dat_iucn_taxonomy$Species,
                                     eco_data$Species,
                                     distrib_sp_iucn))

iucn_dd_total <- Reduce(intersect, list(dat_iucn_dd$Species, dat_iucn_taxonomy_dd$Species,
                                        eco_data_dd$Species,
                                        distrib_sp_iucn_dd))

cas_total <- intersect(fisb_data_cas$Species, new_shp_cas)

total_sp_dd <- c(iucn_dd_total, cas_total) ## 120 DD species with all data

total_sp <- c(iucn_total, iucn_dd_total, cas_total)

## Species to imput
imp_sp <- setdiff(total_sp, shark_trees_clean[[1]]$tip.label) ##104 species

## Imputation
shark_trees_complete <- lapply(shark_trees_clean, function(tree) {
  
  my.info.noranks <- build_info(
    imp_sp, 
    tree = tree,
    find.ranks = FALSE,
    mode = "backbone"
  )
  
  my.input.noranks <- info2input(my.info.noranks, tree)
  
  tree_filled <- rand_tip(
    my.input.noranks, 
    tree,
    rand.type = "random",
    prob = FALSE,
    prune = FALSE,
    forceultrametric = TRUE
  )
  
  return(tree_filled)
})

## Converting the object from a list to a multiPhylo
shark_trees_complete <- do.call(c, shark_trees_complete)

## Salving final phylogenies
write.tree(shark_trees_complete, "Data/Trees/final_trees.tre")

######################################################################################################################################
# Creating a final dataset with intrinsic traits, distribution, and IUCN data and a final set of phylogenies with the same species

# IUCN
common_iucn <- Reduce(intersect, list(dat_iucn$Species, dat_iucn_taxonomy$Species,
                                      eco_data$Species,
                                      distrib_sp_iucn,
                                      shark_trees_complete[[1]]$tip.label)) #979 species

iucn <- dat_iucn[dat_iucn$Species %in% common_iucn, ]
taxons <- dat_iucn_taxonomy[dat_iucn_taxonomy$Species %in% common_iucn, ]
eco <- eco_data[eco_data$Species %in% common_iucn, ]

shp_iucn <- new_shp[new_shp$Species %in% common_iucn, ]
distrib_common <- distrib[distrib$Species %in% common_iucn, ]
if (st_crs(distrib_common) != st_crs(shp_iucn)) {
  shp_iucn  <- st_transform(shp_iucn, st_crs(distrib_common))
}
distrib_final <- bind_rows(distrib_common, shp_iucn)

shark_trees_iucn <- lapply(shark_trees_complete, function(tr) {
  valid_tips <- intersect(common_iucn, tr$tip.label)
  keep.tip(tr, valid_tips)
})
shark_trees_iucn <- do.call(c, shark_trees_iucn)

## Juntar dados
dat <- merge(iucn, taxons, by = "Species")
dat <- merge(dat, eco, by = "Species")

# IUCN DD
common_iucn_dd <- Reduce(intersect, list(dat_iucn_dd$Species, dat_iucn_taxonomy_dd$Species,
                                         eco_data_dd$Species,
                                         distrib_sp_iucn_dd,
                                         shark_trees_complete[[1]]$tip.label))

iucn_dd <- dat_iucn_dd[dat_iucn_dd$Species %in% common_iucn_dd, ]
taxons_dd <- dat_iucn_taxonomy_dd[dat_iucn_taxonomy_dd$Species %in% common_iucn_dd, ]
eco_dd <- eco_data_dd[eco_data$Species %in% common_iucn_dd, ]

shp_iucn_dd <- new_shp[new_shp$Species %in% common_iucn_dd, ]
distrib_common_dd <- distrib_dd[distrib_dd$Species %in% common_iucn_dd, ]

if (st_crs(distrib_common_dd) != st_crs(shp_iucn_dd)) {
  shp_iucn_dd <- st_transform(shp_iucn_dd, st_crs(distrib_common_dd))
}
distrib_final_dd <- bind_rows(distrib_common_dd, shp_iucn_dd)

## Juntar dados DD
dat_dd <- merge(iucn_dd, taxons_dd, by = "Species")
dat_dd <- merge(dat_dd, eco_dd, by = "Species") %>%
  select(Species, Length, `Upper Depth`, `Lower Depth`)

# CAS
cas_common <- intersect(fisb_data_cas$Species, new_shp_cas)

cas_data <- fisb_data_cas[fisb_data_cas$Species %in% cas_common, ] %>%
  as.data.frame()

shp_cas <- new_shp[new_shp$Species %in% cas_common, ]

# Juntar dados DD e CAS
dat_dd_final <- bind_rows(dat_dd, cas_data)

if (st_crs(distrib_final_dd) != st_crs(shp_cas)) {
  shp_cas <- st_transform(shp_cas, st_crs(distrib_final_dd))
}
distrib_dd_final <- bind_rows(distrib_final_dd, shp_cas)

###############################################################################
# Salving

## Salving Chondrichthyes iucn data
write.tree(shark_trees_iucn, "/Applications/Chondrichtyes_article/Data/Trees/chond_iucn.tre")
write.csv(dat, "/Applications/Chondrichtyes_article/Data/Secondary data/Intermediary/dat_iucn_intrinsic.csv")
st_write(distrib_final, "/Applications/Chondrichtyes_article/Data/Secondary data/Intermediary/distrib_final_iucn.gpkg")

## Salving Chondrichthyes iucn DD data
write.csv(dat_dd_final, "/Applications/Chondrichtyes_article/Data/Secondary data/Intermediary/dat_iucn_dd_intrinsic.csv")
st_write(distrib_dd_final, "/Applications/Chondrichtyes_article/Data/Secondary data/Intermediary/distrib_final_iucn_dd.gpkg")
