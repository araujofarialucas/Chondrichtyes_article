library(dplyr)
library(ape)
library(sf)
library(geodata)
library(terra)
library(exactextractr)
library(picante)
library(rfishbase)

# Load Files

distrib_iucn <- st_read("distrib_final_iucn.gpkg")
distrib_iucn <- distrib_iucn[match(unique(distrib_iucn$Species), distrib_iucn$Species), ]

distrib_dd <- st_read("distrib_final_iucn_dd.gpkg")
distrib_dd <- distrib_dd[match(unique(distrib_dd$Species), distrib_dd$Species), ]

dat_iucn <- read.csv("dat_iucn_intrinsic.csv") %>%
  mutate(Species = gsub(" ", "_", Species))

dat_dd <- read.csv("dat_iucn_dd_intrinsic.csv") %>%
  mutate(Species = gsub(" ", "_", Species))

trees <-read.tree("final_trees.tre")

##################################################################################################
# Calculating ED

## Calculating ED for all trees
ed_function <- function(trees) {
  results <- vector("list", length(trees))
  for(i in seq_along(trees)) {
    message("Calculando ED: ", i)
    tree <- trees[[i]]
    results[[i]] <- evol.distinct(tree, type = "fair.proportion")
  }
  return(results)
}

ED <- ed_function(trees)
ED_together <- do.call(rbind, ED)

## Verifying the distribution of EDs for each specie across all trees
ED_distrib <- as.data.frame(tapply(ED_together$w, ED_together$Species,
                                   function(x) shapiro.test(x)$p.value)) %>%
  rename(ED_p = "tapply(ED_together$w, ED_together$Species, function(x) shapiro.test(x)$p.value)")

## Calculating ED mean and median for each specie
means <- aggregate(w ~ Species, data = ED_together, FUN = mean, na.rm = TRUE)
medians <- aggregate(w ~ Species, data = ED_together, FUN = median, na.rm = TRUE)

ED_final <- data.frame(
  Species = rownames(ED_distrib),
  p_value = ED_distrib[, 1],
  Mean = means$w,
  Median = medians$w
)

## Selecting the mean or median ED value for each specie (depending on the distribution)
ED_final$Valor_Final <- ifelse(ED_final$p_value > 0.05,  
                               ED_final$Mean,
                               ED_final$Median)

ED_final <- ED_final %>%
  rename(ED = Valor_Final) %>%
  select(Species, ED)

## Saving the ED values
write.csv(ED_final, "EDs.csv")

##################################################################################################
# Creating environmental dataset

## Loading ocean environmental layers
temp_mean <- bio_oracle(path = tempdir(), var = "Temperature", stat = "Mean", benthic = FALSE, time = "Present")
temp_range <- bio_oracle(path = tempdir(), var = "Temperature", stat = "Range", benthic = FALSE, time = "Present")
sal_mean  <- bio_oracle(path = tempdir(), var = "Salinity", stat = "Mean", benthic = FALSE, time = "Present")
sal_range <- bio_oracle(path = tempdir(), var = "Salinity", stat = "Range", benthic = FALSE, time = "Present")
o2_mean   <- bio_oracle(path = tempdir(), var = "Dissolved.oxygen", stat = "Mean", benthic = FALSE, time = "Present")
o2_range  <- bio_oracle(path = tempdir(), var = "Dissolved.oxygen", stat = "Range", benthic = FALSE, time = "Present")
ph_mean   <- bio_oracle(path = tempdir(), var = "pH", stat = "", benthic = FALSE, time = "Present")

## Alining rasters
temp_range <- project(temp_range, temp_mean)
sal_mean   <- project(sal_mean,   temp_mean)
sal_range  <- project(sal_range,  temp_mean)
o2_mean    <- project(o2_mean,    temp_mean)
o2_range   <- project(o2_range,   temp_mean)
ph_mean    <- project(ph_mean,    temp_mean)

env_rast <- c(temp_mean, temp_range, sal_mean, sal_range, o2_mean, o2_range, ph_mean)

# Giving names to layers
names(env_rast) <- c("temp_mean", "temp_range", "sal_mean", "sal_range", "O2_mean", "O2_range", "ph_mean")

# Extracting environmental values and creating dataset

extract_val_iucn <- exact_extract(env_rast, distrib_iucn, fun = "mean", append_cols = "Species")
dat_env_iucn <- extract_val_iucn %>%
  mutate(Species = factor(Species, levels = distrib_iucn$Species))
dat_env_iucn <- na.omit(dat_env_iucn)


extract_val_dd <- exact_extract(env_rast, distrib_dd, fun = "mean", append_cols = "Species")
dat_env_dd <- extract_val_dd %>%
  mutate(Species = factor(Species, levels = distrib_dd$Species))
dat_env_dd <- na.omit(dat_env_dd)


###############################################################################
#Creating final datasets

## IUCN
common_iucn <- Reduce(intersect, list(dat_env_iucn$Species, dat_iucn$Species, ED_final$Species))

dat_iucn <- dat_iucn[dat_iucn$Species %in% common_iucn, ]
dat_env_iucn <- dat_env_iucn[dat_env_iucn$Species %in% common_iucn, ]
ED_iucn <- ED_final[ED_final$Species %in% common_iucn, ]

dat_analysis_iucn <- merge(dat_iucn, dat_env_iucn, by = "Species")
dat_analysis_iucn <- merge(dat_analysis_iucn, ED_iucn, by = "Species")

### Final adjustments
dat_analysis_iucn$Length <- log(dat_analysis_iucn$Length)
dat_analysis_iucn$IUCN <- ifelse(dat_analysis_iucn$IUCN %in% c("Least Concern", "Near Threatened"), 0, 1)

## DD
common_dd <- Reduce(intersect, list(dat_env_dd$Species, dat_dd$Species, ED_final$Species))

dat_dd <- dat_dd[dat_dd$Species %in% common_dd, ]
dat_env_dd <- dat_env_dd[dat_env_dd$Species %in% common_dd, ]
ED_dd <- ED_final[ED_final$Species %in% common_dd, ]

dat_analysis_dd <- merge(dat_dd, dat_env_dd, by = "Species")
dat_analysis_dd <- merge(dat_analysis_dd, ED_dd, by = "Species")

### Final adjustments
dat_analysis_dd$Length <- log(dat_analysis_dd$Length)

###############################################################################
# Salving

## Salving Chondrichtyes dataset
write.csv(dat_analysis_iucn, "Data/Secondary data/Analysis data/dat_analysis_chond.csv")

## Salving Groups datasets

### Creating sharks and rays datasets
orders_sharks <- c(
  "CARCHARHINIFORMES",
  "HETERODONTIFORMES",
  "HEXANCHIFORMES",
  "LAMNIFORMES",
  "ORECTOLOBIFORMES",
  "PRISTIOPHORIFORMES",
  "SQUALIFORMES",
  "SQUATINIFORMES",
  "ECHINORHINIFORMES"
)

orders_rays <- c(
  "MYLIOBATIFORMES",
  "RAJIFORMES",
  "RHINOPRISTIFORMES",
  "TORPEDINIFORMES"
)


#### IUCN
dat_shark_iucn <- dat_analysis_iucn[dat_analysis_iucn$Order %in% orders_sharks, ]
dat_ray_iucn <- dat_analysis_iucn[dat_analysis_iucn$Order %in% orders_rays, ]

common_shark <- intersect(dat_shark_iucn$Species, trees[[1]]$tip.label)
shark_trees_sub <- lapply(trees, keep.tip, tip = common_shark)
shark_trees_sub <- do.call(c, shark_trees_sub)

common_ray <- intersect(dat_ray_iucn$Species, trees[[1]]$tip.label)
ray_trees_sub <- lapply(trees, keep.tip, tip = common_ray)
ray_trees_sub <- do.call(c, ray_trees_sub)

#### DD

dat_analysis_dd <- dat_analysis_dd %>%
  mutate(Species = gsub("_", " ", Species))

orders_table <- load_taxa() %>% 
  filter(Species %in% dat_analysis_dd$Species) %>% 
  select(Species, Order)

dat_analysis_dd <- merge(dat_analysis_dd, orders_table, by = "Species")

orders_sharks_dd <- c(
  "Carcharhiniformes",
  "Heterodontiformes",
  "Hexanchiformes",
  "Lamniformes",
  "Orectolobiformes",
  "Pristiophoriformes",
  "Squaliformes",
  "Squatiniformes",
  "Echinorhiniformes"
)

orders_rays_dd <- c(
  "Myliobatiformes",
  "Rajiformes",
  "Rhinopristiformes",
  "Torpediniformes"
)

dat_shark_dd <- dat_analysis_dd[dat_analysis_dd$Order %in% orders_sharks_dd, ]
dat_ray_dd <- dat_analysis_dd[dat_analysis_dd$Order %in% orders_rays_dd, ]

### Saving subsets
write.csv(dat_shark_iucn, "Data/Secondary data/Analysis data/shark_data_analysis.csv")
write.csv(dat_ray_iucn, "Data/Secondary data/Analysis data/ray_data_analysis.csv")

write.tree(shark_trees_sub , "Data/Trees/shark_subset.tre")
write.tree(ray_trees_sub, "Data/Trees/ray_subset.tre")

write.csv(dat_shark_dd, "Data/Secondary data/Analysis data/shark_data_predict.csv")
write.csv(dat_ray_dd, "Data/Secondary data/Analysis data/ray_data_predict.csv")

