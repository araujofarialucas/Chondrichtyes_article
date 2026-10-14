# Loading packages
library(ape)
library(caper)
library(phyr)
library(dplyr)

# Loading data
dat <- read.csv("Data/Secondary data/Analysis data/dat_analysis.csv")
chond_trees <- read.tree("Data/Trees/chond.tre")

###############################################################################

# Measuring phylogenetic signal

d <- function(chond_trees, dat) {
  resultados <- vector("list", length(chond_trees))
  for (i in seq_along(chond_trees)) {
    message("Rodando a árvore índice: ", i)
    tr <- chond_trees[[i]]
    tr$node.label <- NULL
    resultados[[i]] <- tryCatch({
      phylo.d(dat, tr, Species, IUCN, permut = 9999)
    }, error = function(e) {
      message("--> ERRO na árvore ", i, ": ", e$message)
      NULL
    })
  }
  return(resultados)
}

d_list <- d(chond_trees, dat)
d_vals <- sapply(d_list, function(x) x$DEstimate)

hist(d_vals,
     probability = TRUE,
     col = "lightblue",
     border = "white",
     main = paste("Histograma - P-values Linha", i), 
     xlab = "Valores", 
     ylab = "Densidade")
curve(dnorm(x, mean = mean(d_vals, na.rm = TRUE), sd = sd(dados, na.rm = TRUE)),
      col = "red",
      lwd = 2,
      add = TRUE)

mean(d_vals, na.rm = TRUE)
median(d_vals, na.rm = TRUE)
quantile(d_vals, probs = c(0.025, 0.975), na.rm = TRUE)

sum(sapply(d_list, function(x) x$Pval1) < 0.05, na.rm = TRUE)
sum(sapply(d_list, function(x) x$Pval0) < 0.05, na.rm = TRUE)


