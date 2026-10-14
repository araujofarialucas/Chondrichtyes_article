# Baixar pacotes
library(ape)
library(caper)
library(dplyr)
library(randomForest)

# Baixar dados
dat_shark <- read.csv("Data/Secondary data/Analysis data/shark_data_analysis.csv")
shark_trees <- read.tree("Data/Trees/shark_subset.tre")
dat_dd <- read.csv("Data/Secondary data/Analysis data/shark_data_predict.csv")

##################################################################################################

# Measuring phylogenetic signal

d_shark <- function(shark_trees, dat_shark) {
  resultados <- vector("list", length(shark_trees))
  for (i in seq_along(shark_trees)) {
    message("Rodando a árvore índice: ", i)
    tr <- shark_trees[[i]]
    tr$node.label <- NULL
    resultados[[i]] <- tryCatch({
      phylo.d(dat_shark, tr, Species, IUCN, permut = 9999)
    }, error = function(e) {
      message("--> ERRO na árvore ", i, ": ", e$message)
      NULL
    })
  }
  return(resultados)
}

d_list_shark <- d_shark(shark_trees, dat_shark)
d_vals_shark <- sapply(d_list_shark, function(x) x$DEstimate)

hist(d_vals_shark,
     probability = TRUE,
     col = "lightblue",
     border = "white",
     main = "Histograma - d_vals_shark", 
     xlab = "Valores", 
     ylab = "Densidade")
curve(dnorm(x, mean = mean(d_vals_shark, na.rm = TRUE), sd = sd(d_vals_shark, na.rm = TRUE)),
      col = "red",
      lwd = 2,
      add = TRUE)

mean(d_vals_shark, na.rm = TRUE)
median(d_vals_shark, na.rm = TRUE)
quantile(d_vals_shark, probs = c(0.025, 0.975), na.rm = TRUE)

sum(sapply(d_list_shark, function(x) x$Pval1) < 0.05, na.rm = TRUE)
sum(sapply(d_list_shark, function(x) x$Pval0) < 0.05, na.rm = TRUE)

# Searching for most threatened orders
table(dat_shark$Order, dat_shark$IUCN)

pct_order <- tapply(dat_shark$IUCN == "1", dat_shark$Order, mean) * 100
pct_order[pct_order > 75]
pct_order[pct_order < 25]

######################################################################################

# Regression

## Looking for colinearity
dat_shark_numeric <- dat_shark[, c("Length", "Upper.Depth", "Lower.Depth", "ED",
                                   "mean.temp_mean","mean.temp_range", "mean.sal_mean",
                                   "mean.sal_range", "mean.O2_mean", "mean.O2_range", "mean.ph_mean")]

matriz_cor <- cor(dat_shark_numeric)

## Runing analysis
rownames(dat_shark) <- dat_shark$Species

set.seed(123)
BPGLM_shark <- function(shark_trees, dat_shark){
  results <- vector("list", length(shark_trees))
  for(i in seq_along(shark_trees)){
    message("rodando árvore: ", i)
    results[[i]] <- binaryPGLMM(IUCN ~ Length + Upper.Depth + Lower.Depth +
                                  mean.temp_mean + mean.temp_range +
                                  mean.sal_mean + mean.ph_mean + ED,
                                data = dat_shark, 
                                phy = shark_trees[[i]])
  }
  return(results)
}

regression_s <- BPGLM_shark(shark_trees, dat_shark)
p_vals_shark <- as.data.frame(sapply(regression_s, function(x) x$B.pvalue))
est_vals_shark <- as.data.frame(sapply(regression_s, function(x) x$B))
s2_vals_shark <- as.data.frame(sapply(regression_s, function(x) x$s2))

### P-values
for (i in 1:nrow(p_vals_shark)) {
  dados <- as.numeric(unlist(p_vals_shark[i, ]))
  hist(dados, 
       probability = TRUE, 
       col = "lightblue", 
       border = "white",
       main = paste("Histograma - P-values Linha", i), 
       xlab = "Valores", 
       ylab = "Densidade")
  curve(dnorm(x, mean = mean(dados, na.rm = TRUE), sd = sd(dados, na.rm = TRUE)), 
        col = "red", 
        lwd = 2, 
        add = TRUE)
  if (i < nrow(p_vals_shark)) {
    readline(prompt = "Pressione [Enter] para ver o próximo gráfico...")
  }
}

which_signif <- function(df) {
  results <- vector("list", nrow(df))
  for(i in seq_len(nrow(df))) {
    linha_dados <- as.numeric(unlist(df[i, ]))
    results[[i]] <- sum(linha_dados < 0.05, na.rm = TRUE)
  }
  return(results)
}
which_signif(p_vals_shark)

mean_p_s <- function(df) {
  results <- vector("list", nrow(df))
  for(i in seq_len(nrow(df))) {
    linha_dados <- as.numeric(unlist(df[i, ]))
    results[[i]] <- sum(linha_dados, na.rm = TRUE) / ncol(df)
  }
  return(results)
}
mean_p_s(p_vals_shark)

p_sd <- function(df) {
  results <- vector("list", nrow(df))
  for(i in seq_len(nrow(df))) {
    dados <- as.numeric(unlist(df[i, ]))
    results[[i]] <- quantile(dados, probs = c(0.025, 0.975), na.rm = TRUE)
  }
  return(results)
}
p_sd(p_vals_shark)

### Estimate values
for (i in 1:nrow(est_vals_shark)) {
  dados <- as.numeric(unlist(est_vals_shark[i, ]))
  hist(dados, 
       probability = TRUE, 
       col = "lightblue", 
       border = "white",
       main = paste("Histograma - Estimativas Linha", i), 
       xlab = "Valores", 
       ylab = "Densidade")
  curve(dnorm(x, mean = mean(dados, na.rm = TRUE), sd = sd(dados, na.rm = TRUE)), 
        col = "red", 
        lwd = 2, 
        add = TRUE)
  if (i < nrow(est_vals_shark)) {
    readline(prompt = "Pressione [Enter] para ver o próximo gráfico...")
  }
}

mean_est <- function(df) {
  results <- vector("list", nrow(df))
  for(i in seq_len(nrow(df))) {
    linha_dados <- as.numeric(unlist(df[i, ]))
    results[[i]] <- sum(linha_dados, na.rm = TRUE) / ncol(df)
  }
  return(results)
}
mean_est(est_vals_shark)

est_sd <- function(df) {
  results <- vector("list", nrow(df))
  for(i in seq_len(nrow(df))) {
    dados <- as.numeric(unlist(df[i, ]))
    results[[i]] <- quantile(dados, probs = c(0.025, 0.975), na.rm = TRUE)
  }
  return(results)
}
est_sd(est_vals_shark)

### s2 values
s2_vec_s <- unlist(sapply(regression_s, function(x) x$s2))
hist(s2_vec_s, probability = TRUE, main = "Histograma de s2", xlab = "s2", col = "lightblue", border = "white")
curve(dnorm(x, mean = mean(s2_vec_s), sd = sd(s2_vec_s)), 
      col = "red", lwd = 2, add = TRUE)

mean_s2_s <- lapply(s2_vals_shark, mean)
median_s2_s <- lapply(s2_vals_shark, median)
CI_s2_s <- lapply(s2_vals_shark, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))

############################################################################################################

# Random forest

## Selecting best parameter values
dat_shark$IUCN <- as.factor(dat_shark$IUCN)
table(dat_shark$IUCN)

for (ns in c(1, 5, 10, 20, 30, 40, 50)) {
  for (p in seq(1.0, 1.8, by = 0.1)) {
    mod <- randomForest(IUCN ~ Length + Upper.Depth + Lower.Depth +
                          mean.temp_mean + mean.temp_range +
                          mean.sal_mean + mean.ph_mean + ED,
                        data = dat_shark, mtry = 2, nodesize = ns, ntree = 500,
                        classwt = c("0" = 1.0, "1" = p))
    a0 <- 1 - mod$confusion["0", "class.error"]
    a1 <- 1 - mod$confusion["1", "class.error"]
    acc_inv <- (a0 * (165 / 464)) + (a1 * (299 / 464))
    cat("Node:", ns, "| Peso:", p, "| Acerto Ponderado:", round(acc_inv * 100, 1), "%\n")
  }
}

## Runing model
set.seed(123)
rf_modelo_s <- randomForest(IUCN ~ Length + Upper.Depth + Lower.Depth +
                            mean.temp_mean + mean.temp_range +
                            mean.sal_mean + mean.ph_mean + ED,
                          data = dat_shark,
                          mtry = 2, replace = TRUE, nodesize = 40,
                          ntree = 500,
                          importance = TRUE, do.trace = TRUE,
                          classwt = c("0" = 1.0, "1" = 1.7))

## Looking at traits importance and model prediction accuracy
importance(rf_modelo_s)
rf_modelo_s

## Predicting threat risk for not assessed species
set.seed(123)
dat_dd$pred_IUCN <- predict(rf_modelo_s, newdata = dat_dd, type = "response")

df_predicoes <- data.frame(
  Species = dat_dd$Species,
  Pred_IUCN = dat_dd$pred_IUCN
)

table(df_predicoes$Pred_IUCN)

write.csv(df_predicoes, "Data/Secondary data/Results/Prediction_sharks.csv")
############################################################################################################



