# ------------------------ Calculando o Range ----------------------------#
library(sf)
library(rgbif)
library(tidyverse)
library(letsR)
library(dbscan)

dados1 <- read.csv("Data/Secondary data/Manual inspection/Species without shp/cas_species.csv")
dados2 <- read.csv("Data/Secondary data/Manual inspection/Species without shp/dd_without_shp.csv")
dados3 <- read.csv("Data/Secondary data/Manual inspection/Species without shp/iucn_without_shp.csv")


carentes <- rbind.data.frame(dados1, dados2, dados3)
carentes[1] <- 1:length(carentes$X)


get_m_occs <- function(vetor){
  a <- vector(mode = 'list', length = length(vetor))
  
  for(i in 1:length(vetor)){
    
    print(i)
    
    b <- occ_search(scientificName = vetor[i])
    b <- b$data
    if("decimalLongitude" %in% colnames(b) == FALSE | 'species' %in% colnames(b) == FALSE){next}
    else{b <- b%>%
      dplyr::select(species, decimalLongitude, decimalLatitude)}
    
    a[[i]] <- b
  }
  
  a <- do.call('rbind.data.frame', a)
  
  a <- a%>%
    dplyr::filter(!is.na(decimalLongitude) | !is.na(decimalLatitude))
  names(a) <- c('species', 'longitude', 'latitude')
  
  return(a)
  
}



especies <- carentes$x

tubaroes <- get_m_occs(especies)

names(tubaroes) <-c("scientificName", "decimalLongitude", "decimalLatitude")


tubaroes <- tubaroes%>%
  filter(!is.na(decimalLongitude) | !is.na(decimalLatitude))

view(tubaroes)


#Necessário pq a função só aceita assim


#Obtendo shp p/ intersecto
agua <- st_read(
  "Data/Primary data/World Bank Official Boundaries - Ocean Mask"
)

#link p/ shp:
#https://datacatalog.worldbank.org/search/dataset/0038272/world-bank-official-boundaries

#Esse shapefile é útil mas não tem uma qualidade muito boa... Se conseguir
#achar um com uma qualidade melhor, é o ideal



#Aplicando o lets_range



ranges <- lets.range(occs = tubaroes, crop_shp = agua, buf = 500000,
                   ratio = 1,
                     intersection = TRUE, eps_m = 2000000)

st_write(ranges, 'Shapefiles/shark_ranges.shp')

sharks <- st_read('Shapefiles/Ranges/shark_ranges.shp')
