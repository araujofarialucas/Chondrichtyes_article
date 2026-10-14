# Loading packages
library(ape)
library(dplyr)
library(stringr)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(tidyverse)
library(ggforce)
library(gridExtra)

# Loading data
dat_analysis <- read.csv("Data/Secondary data/Final data/dat_analysis.csv") %>%
  mutate(Ordem = str_to_sentence(Ordem)) %>%
  dplyr::select(Species, Ordem, IUCN) %>%
  rename(Order = Ordem)
rownames(dat_analise) <- dat_analise$Species

chond_tree <- read.tree("Data/Trees/chond.tre")

###############################################################################

# Ploting tree

orders_sharks <- c(
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

orders_rays <- c(
  "Myliobatiformes",
  "Rajiformes",
  "Rhinopristiformes",
  "Torpediniformes"
)

dat_analysis$Group <- ifelse(dat_analysis$Order %in% ordens_sharks, "Shark",
                    ifelse(dat_analysis$Order %in% ordens_rays, "Ray", "Chimaera"))

dat_analysis$Risk <- ifelse(dat_analysis$IUCN == 0, "Non-threatened", "Threatened")

collors_manual_orders <- c(
  "Carcharhiniformes"  = "#8E1EA2",
  "Chimaeriformes"     = "#377eb8",
  "Heterodontiformes"  = "#4daf4a", 
  "Hexanchiformes"     = "#AEC4D4", 
  "Lamniformes"        = "#ff7f00", 
  "Myliobatiformes"    = "#ffff33", 
  "Orectolobiformes"   = "#a65628",
  "Pristiophoriformes" = "#4D2A00", 
  "Rajiformes"         = "#999999",
  "Rhinopristiformes"  = "#66c2a5", 
  "Squaliformes"       = "#fc8d62", 
  "Squatiniformes"     = "#8da0cb",
  "Torpediniformes"    = "#e78ac3" 
)

collors_manual_threat <- c(
  "Non-threatened" = "#228B22", "Threatened" = "#D81E05")


p <- ggtree(tree_unica, layout = "circular", linewidth = 0.2, open.angle = 10)

# Creating first ring (threat status)
p_ring1 <- gheatmap(
  p,
  dat_analise %>% dplyr::select(Risk), 
  offset = 0.005,   
  width = 0.05,    
  colnames = FALSE,
  color = NA
) +
  scale_fill_manual(
    values = cores_manual_threat,
    name = NULL) +
  theme(legend.position = "none")

# Creating a new scale
p_ring1 <- p_ring1 + new_scale_fill()

# Creating second ring (orders)
p_ring2 <- gheatmap(
  p_ring1,
  dat_analise %>% dplyr::select(Order), 
  offset = 15,
  width = 0.05,    
  colnames = FALSE,
  color = NA
) + scale_fill_manual(
  values = cores_manual_ordens,
  name = "Orders"
  ) + 
  theme(legend.position = "none")

ggsave("grafico_perfeito.png", plot = p_ring2, dpi = 600, width = 10, height = 10, units = "in", bg = "transparent")

###############################################################################
# Regressions

dat_s <- read.csv("Data/Secondary data/Final data/shark_data.csv")
dat_r <- read.csv("Data/Secondary data/Final data/ray_data.csv")

dat_s$Grupo <- "Sharks"
dat_r$Grupo <- "Rays"
dat_total <- bind_rows(dat_s, dat_r)

dat_total$IUCN_graph <- ifelse(dat_total$IUCN == 0, "Non-Threatened", "Threatened")

dat_combined_plot <- dat_total %>%
  pivot_longer(
    cols = c("Length", "mean.sal_mean", "mean.ph_mean", "Upper.Depth",
             "Lower.Depth", "mean.temp_mean", "mean.temp_range", "ED"),
    names_to = "Variavel", 
    values_to = "Valor"
  ) %>%
  dplyr::select(Species, Grupo, IUCN_graph, Variavel, Valor) %>%
  mutate(
    Panel = case_when(
      Variavel == "Length" ~ "Body Length",
      Variavel == "Upper.Depth" ~ "Upper Depth",
      Variavel == "Lower.Depth" ~ "Lower Depth",
      Variavel == "mean.ph_mean" ~ "Mean pH",
      Variavel == "mean.sal_mean" ~ "Mean Salinity",
      Variavel == "mean.temp_mean" ~ "Mean Temperature",
      Variavel == "mean.temp_range" ~ "Temperature Range",
      Variavel == "ED" ~ "ED",
      TRUE ~ Variavel
    )
  )

sig_data <- data.frame(
  Panel = rep(c("Body Length", "Upper Depth", "Lower Depth", "Mean pH", 
                "Mean Salinity", "Mean Temperature", "Temperature Range", "ED"), 
              each = 2),
  Grupo = rep(c("Sharks", "Rays"), times = 8),
  Sig = c(
    "**", "**",  # Body Length (Sharks, Rays)
    "", "*",      # Upper Depth
    "**", "**",      # Lower Depth
    "**", "**",  # Mean pH
    "**", "**",  # Mean Salinity
    "", "",      # Mean Temperature
    "", "",     # Temperature Range
    "", ""       # ED
  ),
  y_pos = rep(c(
    7.5,    # Body Length
    600,    # Upper Depth
    2500,   # Lower Depth
    8.4,    # Mean pH
    38,     # Mean Salinity
    30,     # Mean Temperature
    22,     # Temperature Range
    60      # ED
  ), each = 2)
)

dat_lines <- dat_combined_plot %>%
  group_by(Panel, Grupo, IUCN_graph) %>%
  summarise(med = median(Valor, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = IUCN_graph, values_from = med) %>%
  mutate(
    x_num = as.numeric(factor(Grupo)),
    x_start = x_num - 0.175, 
    x_end   = x_num + 0.175  
  )

p_single <- ggplot(dat_combined_plot, aes(x = Grupo, y = Valor, fill = IUCN_graph, color = IUCN_graph)) +
  geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.3, 
               position = position_dodge(0.7)) +
  
  geom_segment(
    data = dat_lines,
    aes(x = x_start, y = `Non-Threatened`, xend = x_end, yend = Threatened),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.8
  ) +
  
  scale_fill_manual(values = c(
    "Non-Threatened" = "#228B22", 
    "Threatened"     = "#D81E05"  
  )) +
  scale_color_manual(values = c(
    "Non-Threatened" = "#228B22", 
    "Threatened"     = "#D81E05"  
  )) +
  geom_text(
    data = sig_data, 
    aes(x = Grupo, y = y_pos, label = Sig), 
    inherit.aes = FALSE, 
    position = position_dodge(0.7),
    size = 8, fontface = "bold", color = "black"
  ) +
  
  facet_wrap(~ Panel, scales = "free", ncol = 2) +
  labs(title = NULL, x = NULL, y = NULL) + 
  theme_classic(base_size = 11) + 
  theme(
    axis.text.y = element_text(size = 21, color = "black"),
    axis.text.x = element_text(size = 21, color = "black", face = "bold"), 
    strip.text = element_text(size = 21, color = "black", face = "bold", angle = 0, hjust = 0),
    strip.background = element_blank(),
    legend.position = "none",
    panel.spacing.x = unit(2.5, "lines"), 
    panel.spacing.y = unit(2, "lines"),
    plot.margin = margin(t = 15, r = 20, b = 15, l = 15, unit = "pt")
  )

print(p_single)

ggsave("Figures/sharks_rays_single_figure.jpeg", plot = p_single, width = 16, height = 15, dpi = 300, limitsize = FALSE)

