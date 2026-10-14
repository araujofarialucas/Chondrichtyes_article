# Loading packages
library(ape)
library(caper)
library(dplyr)
library(randomForest)

# Loading data
dat_ray <- read.csv("Data/Secondary data/Analysis data/ray_data_analysis.csv")
ray_trees <- read.tree("Data/Trees/ray_subset.tre")
dat_dd <- read.csv("Data/Secondary data/Analysis data/ray_data_predict.csv")

####################################################################################################

# Measuring phylogenetic signal

d_ray <- function(ray_trees, dat_ray) {
  resultados <- vector("list", length(ray_trees))
  for (i in seq_along(ray_trees)) {
    message("Rodando a árvore índice: ", i)
    tr <- ray_trees[[i]]
    tr$node.label <- NULL
    resultados[[i]] <- tryCatch({
      phylo.d(dat_ray, tr, Species, IUCN, permut = 9999)
    }, error = function(e) {
      message("--> ERRO na árvore ", i, ": ", e$message)
      NULL
    })
  }
  return(resultados)
}

d_list_ray <- d_ray(ray_trees, dat_ray)
d_vals_ray <- sapply(d_list_ray, function(x) x$DEstimate)

hist(d_vals_ray,
     probability = TRUE,
     col = "lightblue",
     border = "white",
     main = paste("Histograma - P-values Linha", i), 
     xlab = "Valores", 
     ylab = "Densidade")
curve(dnorm(x, mean = mean(d_vals_ray, na.rm = TRUE), sd = sd(dados, na.rm = TRUE)),
      col = "red",
      lwd = 2,
      add = TRUE)

mean(d_vals_ray , na.rm = TRUE)
median(d_vals_ray , na.rm = TRUE)
quantile(d_vals_ray, probs = c(0.025, 0.975), na.rm = TRUE)

sum(sapply(d_list_ray, function(x) x$Pval1) < 0.05, na.rm = TRUE)
sum(sapply(d_list_ray, function(x) x$Pval0) < 0.05, na.rm = TRUE)

# Searching for most threatened orders
table(dat_ray$Order, dat_ray$IUCN)

pct_order <- tapply(dat_ray$IUCN == 1, dat_ray$Order, mean, na.rm = TRUE) * 100
pct_order[pct_order > 75]
pct_order[pct_order < 25]

######################################################################################

# Regression

## Looking for colinearity
dat_ray_numeric <- dat_ray[, c("Length", "Upper.Depth", "Lower.Depth", "ED",
                               "mean.temp_mean","mean.temp_range", "mean.sal_mean",
                               "mean.sal_range", "mean.O2_mean", "mean.O2_range", "mean.ph_mean")]

matriz_cor <- cor(dat_ray_numeric)

## Runing analysis
rownames(dat_ray) <- dat_ray$Species

BPGLM_ray <- function(ray_trees, dat_ray){
  results <- vector("list", length(ray_trees))
  for(i in seq_along(ray_trees)){
    message("rodando árvore: ", i)
    results[[i]] <- binaryPGLMM(IUCN ~ Length + Upper.Depth + Lower.Depth +
                                  mean.temp_mean + mean.temp_range +
                                  mean.sal_mean + mean.ph_mean + ED,
                                data = dat_ray, 
                                phy = ray_trees[[i]])
  }
  return(results)
}

regression_r <- BPGLM_ray(ray_trees, dat_ray)
p_vals_ray <- as.data.frame(sapply(regression_r, function(x) x$B.pvalue))
est_vals_ray <- as.data.frame(sapply(regression_r, function(x) x$B))
s2_vals_ray <- as.data.frame(sapply(regression_r, function(x) x$s2))

### P-values
for (i in 1:nrow(p_vals_ray)) {
  dados <- as.numeric(unlist(p_vals_ray[i, ]))
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
  if (i < nrow(p_vals_ray)) {
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
which_signif(p_vals_ray)

mean_p_r <- function(df) {
  results <- vector("list", nrow(df))
  for(i in seq_len(nrow(df))) {
    linha_dados <- as.numeric(unlist(df[i, ]))
    results[[i]] <- sum(linha_dados, na.rm = TRUE) / ncol(df)
  }
  return(results)
}
mean_p_r(p_vals_ray)

p_sd_r <- function(df) {
  results <- vector("list", nrow(df))
  for(i in seq_len(nrow(df))) {
    dados <- as.numeric(unlist(df[i, ]))
    results[[i]] <- quantile(dados, probs = c(0.025, 0.975), na.rm = TRUE)
  }
  return(results)
}
p_sd_r(p_vals_ray)

### Estimate values
for (i in 1:nrow(est_vals_ray)) {
  dados <- as.numeric(unlist(est_vals_ray[i, ]))
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
  if (i < nrow(est_vals_ray)) {
    readline(prompt = "Pressione [Enter] para ver o próximo gráfico...")
  }
}

mean_est_r <- function(df) {
  results <- vector("list", nrow(df))
  for(i in seq_len(nrow(df))) {
    linha_dados <- as.numeric(unlist(df[i, ]))
    results[[i]] <- sum(linha_dados, na.rm = TRUE) / ncol(df)
  }
  return(results)
}
mean_est_r(est_vals_ray)

est_sd_r <- function(df) {
  results <- vector("list", nrow(df))
  for(i in seq_len(nrow(df))) {
    dados <- as.numeric(unlist(df[i, ]))
    results[[i]] <- quantile(dados, probs = c(0.025, 0.975), na.rm = TRUE)
  }
  return(results)
}
est_sd_r(est_vals_ray)

### s2 values
s2_vec_r <- unlist(sapply(regression_r, function(x) x$s2))
hist(s2_vec_r, probability = TRUE, main = "Histograma de s2", xlab = "s2", col = "lightblue", border = "white")
curve(dnorm(x, mean = mean(s2_vec_r), sd = sd(s2_vec_r)), 
      col = "red", lwd = 2, add = TRUE)

mean_s2_r <- lapply(s2_vals_ray, mean)
median_s2_r <- lapply(s2_vals_ray, median)
CI_s2_r <- lapply(s2_vals_ray, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))

############################################################################################################

# Random forest

## Selecting best parameter values
dat_ray$IUCN <- as.factor(dat_ray$IUCN)

set.seed(123)
for (ns in c(1, 5, 10, 20, 30, 40, 50)) {
  for (p in seq(1.0, 1.7, by = 0.1)) {
    mod <- randomForest(IUCN ~ Length + Upper.Depth + Lower.Depth +
                          mean.temp_mean + mean.temp_range +
                          mean.sal_mean + mean.ph_mean + ED,
                        data = dat_ray, mtry = 2, nodesize = ns, ntree = 500,
                        classwt = c("0" = 1.0, "1" = p))
    a0 <- 1 - mod$confusion["0", "class.error"]
    a1 <- 1 - mod$confusion["1", "class.error"]
    acc_inv <- (a0 * (175 / 473)) + (a1 * (298 / 473))
    cat("Node:", ns, "| Peso:", p, "| Acerto Ponderado:", round(acc_inv * 100, 1), "%\n")
  }
}

## Runing model
set.seed(123)
rf_modelo_r <- randomForest(IUCN ~ Length + Upper.Depth + Lower.Depth +
                            mean.temp_mean + mean.temp_range +
                            mean.sal_mean + mean.ph_mean + ED,
                          data = dat_ray,
                          mtry = 2, replace = TRUE, nodesize = 30,
                          ntree = 500,
                          importance = TRUE, do.trace = TRUE,
                          classwt = c("0" = 1.0, "1" = 1.1))

## Looking at traits importance and model prediction accuracy
importance(rf_modelo_r)
rf_modelo_r

## Predicting threat risk for not assessed species
set.seed(123)
dat_dd$pred_IUCN <- predict(rf_modelo_r, newdata = dat_dd, type = "response")

df_predicoes <- data.frame(
  Species = dat_dd$Species,
  Pred_IUCN = dat_dd$pred_IUCN
)

table(df_predicoes$Pred_IUCN)

write.csv(df_predicoes, "Data/Secondary data/Results/Prediction_rays.csv")

############################################################################################################

