# ==============================================================================
# ===== (2026) Drosophilidae 16S Virus Infection: Figure - NMDS ================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script runs and plots a NMDS on microbiome ASV relative
# abundances across Drosophilidae host species. Underneath these plots,
# estimates of phylogenetic heritability and repeatability are calculated from
# MCMCglmm models and plotted as mean and HPD interval, with an arbitrary MCMCp
# measure of significance against a threshold of 0.05 used as a test for non-
# zero estimates.

# Due to contamination in samples of D. takahashii and Z. davidi, apparent in
# RNAseq runs of the same samples used for 16S analysis, these species have been
# removed.

# ----- 0.2. Dependencies ------------------------------------------------------

library(tidyverse)
library(MCMCglmm)
library(phyloseq)
library(microbiome)
library(here)
library(ape)
library(patchwork)
library(vegan)
library(zCompositions)
library(compositions)

library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# ----- 0.3. Load Data ---------------------------------------------------------

load(here("data", "data_phyloseq.RData"))

data_metadata <- read_csv(here("data", "data_metadata.csv"))

data_fly <- read_csv(here("data", "data_drosophilidae.csv"))

tree_drosophilidae <- read.tree(here("data", "tree_hosts.nwk"))

# ------------------------------------------------------------------------------
# ----- 1. Data Wrangling ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1 Aggregate by Genus -------------------------------------------------

#data_abundances <- aggregate_taxa(phyloseq_cutoff, level="Genus")

# ----- 1.2 Remove non-experimental samples ------------------------------------

data_abundances <- subset_samples(phyloseq_cutoff, !(ID %in% c(as.character(307:324), "Negative")))
data_abundances <- subset_samples(data_abundances, !(Spc.Name %in% c("D. takahashii", "Z. davidi")))

# ----- 1.3. CLR Transform -----------------------------------------------------

data_clr <- as(otu_table(data_abundances), "matrix") %>%
  cmultRepl(method = "CZM", output = "p-counts", z.warning = 1, z.delete = F)

data_clr <- clr(data_clr) %>% t()

# ----- 1.4. Wrangle to NMDS matrix ---------------------------------------------

rownames(data_clr) <- sample_names(data_abundances)

data_NMDS <- as.data.frame(data_clr)

# ------------------------------------------------------------------------------
# ----- 2. NMDS Analysis & Plots ------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1. NMDS ---------------------------------------------------------------

# sample.int(100000, 1) returned 35621
set.seed(35621)

NMDS_results <- metaMDS(vegdist(as.matrix(data_NMDS), method = "euclidean"), k = 4, trymax = 100)

# ----- 2.2. Join to Metadata --------------------------------------------------

NMDS_results <- as.data.frame(NMDS_results$points)
NMDS_results$sample <- rownames(NMDS_results)

data_abundances <- as.data.frame(as.matrix(sample_data(data_abundances)))

data_abundances <- left_join(data_abundances, NMDS_results)

data_abundances$Condition <- factor(data_abundances$Condition, levels = c("None", "Ringers", "DCV"))

colnames(data_abundances)[21:24] <- c("NMDS1", "NMDS2", "NMDS3", "NMDS4")

# ----- 2.4. Define Theme ------------------------------------------------------

theme_NMDS <- theme(aspect.ratio = 1,
                    text = element_text(size = 8, color = "#2e3440"),
                    panel.grid = element_blank(),
                    legend.position = "none",
                    strip.background = element_rect(fill = "#e5e9f0"))

# ----- 2.5. NMDS1/NMDS2 Plot ------------------------------------------------------

p_NMDS12 <- ggplot(data_abundances) +
  geom_point(aes(x = NMDS1, y = NMDS2, color = Condition), size = 1) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Condition), size = 0.75) +
  scale_x_continuous(name = "NMDS 1") +
  scale_y_continuous(name = "NMDS 2") +
  scale_color_manual(labels = c("Uninjured", "Saline Injection", "DCV Infection"), values = c("#3E356B", "#357BA2", "#49C1AD")) +
  theme_bw() +
  theme_NMDS +
  theme(legend.position = "none") +
  ggtitle("A.")

p_NMDS34 <- ggplot(data_abundances) +
  geom_point(aes(x = NMDS3, y = NMDS4, color = Condition), size = 1) +
  stat_ellipse(aes(x = NMDS3, y = NMDS4, color = Condition), size = 0.75) +
  scale_x_continuous(name = "NMDS 3") +
  scale_y_continuous(name = "NMDS 4") +
  scale_color_manual(labels = c("Uninjured", "Saline Injection", "DCV Infection"), values = c("#3E356B", "#357BA2", "#49C1AD")) +
  theme_bw() +
  theme_NMDS +
  ggtitle("")

pNMDS <- p_NMDS12 | p_NMDS34

# ------------------------------------------------------------------------------
# ----- 3. NMDS Models ----------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 3.1 Thresholded MCMCp function for Heritability/Repeatability ----------

MCMCp <- function(posterior, threshold = 0.05){
  
  p_above <- mean(posterior > threshold)
  p_below <- mean(posterior <= threshold)
  
  p_value <- 2 * min(p_above, p_below)
  
  return(p_value)
  
}

# ----- 3.2 MCMC summary function for Heritability/Repeatability ---------------

summarise_MCMC <- function(quant, metric, data_list) {
  conditions <- c("Uninjured", "Saline Injection", "DCV Infection")
  map_dfr(seq_along(conditions), function(i) {
    data.frame(
      Metric = metric,
      Condition = conditions[i],
      Quant = quant,
      Mean = mean(data_list[[i]]),
      HPD_Low = HPDinterval(data_list[[i]])[1],
      HPD_High = HPDinterval(data_list[[i]])[2],
      p_value = MCMCp(data_list[[i]])
    )
  })
}

# ----- 3.3 Heritability Prior -------------------------------------------------

prior_her <- list(G = list(
  G1 = list(V = diag(3), nu = 3, alpha.mu = rep(0,3), alpha.V = diag(3) * 1000),
  G2 = list(V = diag(3), nu = 3, alpha.mu = rep(0,3), alpha.V = diag(3) * 1000),
  G3 = list(V = diag(1), nu = 0.002),
  G3 = list(V = diag(1), nu = 0.002)),
  R = list(V = diag(3), nu = 0.002))

# ----- 3.4 Repeatability Prior ------------------------------------------------

prior_rep <- list(G = list(
  G1 = list(V = diag(3), nu = 3, alpha.mu = rep(0,3), alpha.V = diag(3) * 1000),
  G2 = list(V = diag(1), nu = 0.002),
  G3 = list(V = diag(1), nu = 0.002)),
  R = list(V = diag(3), nu = 0.002))

# ----- 3.5 Iteration Multiplier -----------------------------------------------

itt <- 5 # Set to 10 for demonstration, 100 for publication

# ----- 3.6 Add columns for phylogenetic MCMCglmms -----------------------------

data_abundances$animal <- data_abundances$Spc.Name %>% str_remove_all(" ")

data_abundances <- as.data.frame(data_abundances)

data_abundances <- left_join(data_abundances, data_fly)

# ----- 3.7 Drop tips for species not in dataset -------------------------------

tree_drosophilidae <- drop.tip(tree_drosophilidae, setdiff(tree_drosophilidae$tip.label, data_abundances$animal))

plot(tree_drosophilidae)

# ----- 3.8 NMDS1 Heritability Model ---------------------------------------------

if (!file.exists(here("models", "ASV_NMDS1_Heritability.RData"))) {
  
  model_NMDS1_her <- MCMCglmm(NMDS1 ~ Condition,
                              random = ~ idh(Condition):animal + idh(Condition):Spc.Name + Diet + ID,
                              rcov = ~idh(Condition):units,
                              pedigree = tree_drosophilidae,
                              prior = prior_her,
                              data = data_abundances,
                              nitt = 130000*itt,
                              thin = 20*itt,
                              burnin = 30000*itt,
                              pr = TRUE)
  
  save(model_NMDS1_her, file = here("models", "ASV_NMDS1_Heritability.RData"))
  
} else {load(here("models", "ASV_NMDS1_Heritability.RData"))}

NMDS1_uninjured_her <- model_NMDS1_her$VCV[,1] / (model_NMDS1_her$VCV[,1] + model_NMDS1_her$VCV[,4])
NMDS1_saline_her <- model_NMDS1_her$VCV[,2] / (model_NMDS1_her$VCV[,2] + model_NMDS1_her$VCV[,5])
NMDS1_dcv_her <- model_NMDS1_her$VCV[,3] / (model_NMDS1_her$VCV[,3] + model_NMDS1_her$VCV[,6])

# ----- 3.9 NMDS1 Repeatability Model --------------------------------------------

if (!file.exists(here("models", "ASV_NMDS1_Repeatability.RData"))) {
  
  model_NMDS1_rep <- MCMCglmm(NMDS1 ~ Condition,
                              random = ~ idh(Condition):animal+ Diet + ID,
                              rcov = ~idh(Condition):units,
                              pedigree = tree_drosophilidae,
                              prior = prior_rep,
                              data = data_abundances,
                              nitt = 130000*itt,
                              thin = 20*itt,
                              burnin = 30000*itt,
                              pr = TRUE)
  
  save(model_NMDS1_rep, file = here("models", "ASV_NMDS1_Repeatability.RData"))
  
} else {load(here("models", "ASV_NMDS1_Repeatability.RData"))}

NMDS1_uninjured_rep <- model_NMDS1_rep$VCV[,1] / (model_NMDS1_rep$VCV[,1] + model_NMDS1_rep$VCV[,6])
NMDS1_saline_rep <- model_NMDS1_rep$VCV[,2] / (model_NMDS1_rep$VCV[,2] + model_NMDS1_rep$VCV[,7])
NMDS1_dcv_rep <- model_NMDS1_rep$VCV[,3] / (model_NMDS1_rep$VCV[,3] + model_NMDS1_rep$VCV[,8])

# ----- 3.10 NMDS2 Heritability Model --------------------------------------------

if (!file.exists(here("models", "ASV_NMDS2_Heritability.RData"))) {
  
  model_NMDS2_her <- MCMCglmm(NMDS2 ~ Condition,
                              random = ~ idh(Condition):animal + idh(Condition):Spc.Name + Diet + ID,
                              rcov = ~idh(Condition):units,
                              pedigree = tree_drosophilidae,
                              prior = prior_her,
                              data = data_abundances,
                              nitt = 130000*itt,
                              thin = 20*itt,
                              burnin = 30000*itt,
                              pr = TRUE)
  
  save(model_NMDS2_her, file = here("models", "ASV_NMDS2_Heritability.RData"))
  
} else {load(here("models", "ASV_NMDS2_Heritability.RData"))}

NMDS2_uninjured_her <- model_NMDS2_her$VCV[,1] / (model_NMDS2_her$VCV[,1] + model_NMDS2_her$VCV[,4])
NMDS2_saline_her <- model_NMDS2_her$VCV[,2] / (model_NMDS2_her$VCV[,2] + model_NMDS2_her$VCV[,5])
NMDS2_dcv_her <- model_NMDS2_her$VCV[,3] / (model_NMDS2_her$VCV[,3] + model_NMDS2_her$VCV[,6])

# ----- 3.11 NMDS2 Repeatability Model -------------------------------------------

if (!file.exists(here("models", "ASV_NMDS2_Repeatability.RData"))) {
  
  model_NMDS2_rep <- MCMCglmm(NMDS2 ~ Condition,
                              random = ~ idh(Condition):animal+ Diet + ID,
                              rcov = ~idh(Condition):units,
                              pedigree = tree_drosophilidae,
                              prior = prior_rep,
                              data = data_abundances,
                              nitt = 130000*itt,
                              thin = 20*itt,
                              burnin = 30000*itt,
                              pr = TRUE)
  
  save(model_NMDS2_rep, file = here("models", "ASV_NMDS2_Repeatability.RData"))
  
} else {load(here("models", "ASV_NMDS2_Repeatability.RData"))}

NMDS2_uninjured_rep <- model_NMDS2_rep$VCV[,1] / (model_NMDS2_rep$VCV[,1] + model_NMDS2_rep$VCV[,6])
NMDS2_saline_rep <- model_NMDS2_rep$VCV[,2] / (model_NMDS2_rep$VCV[,2] + model_NMDS2_rep$VCV[,7])
NMDS2_dcv_rep <- model_NMDS2_rep$VCV[,3] / (model_NMDS2_rep$VCV[,3] + model_NMDS2_rep$VCV[,8])


# ----- 3.12 NMDS3 Heritability Model --------------------------------------------

if (!file.exists(here("models", "ASV_NMDS3_Heritability.RData"))) {
  
  model_NMDS3_her <- MCMCglmm(NMDS3 ~ Condition,
                              random = ~  idh(Condition):animal + idh(Condition):Spc.Name + Diet + ID,
                              rcov = ~idh(Condition):units,
                              pedigree = tree_drosophilidae,
                              prior = prior_her,
                              data = data_abundances,
                              nitt = 130000*itt,
                              thin = 20*itt,
                              burnin = 30000*itt,
                              pr = TRUE)
  
  save(model_NMDS3_her, file = here("models", "ASV_NMDS3_Heritability.RData"))
  
} else {load(here("models", "ASV_NMDS3_Heritability.RData"))}

NMDS3_uninjured_her <- model_NMDS3_her$VCV[,1] / (model_NMDS3_her$VCV[,1] + model_NMDS3_her$VCV[,4])
NMDS3_saline_her <- model_NMDS3_her$VCV[,2] / (model_NMDS3_her$VCV[,2] + model_NMDS3_her$VCV[,5])
NMDS3_dcv_her <- model_NMDS3_her$VCV[,3] / (model_NMDS3_her$VCV[,3] + model_NMDS3_her$VCV[,6])

# ----- 3.13 NMDS3 Repeatability Model -------------------------------------------

if (!file.exists(here("models", "ASV_NMDS3_Repeatability.RData"))) {
  
  model_NMDS3_rep <- MCMCglmm(NMDS3 ~ Condition,
                              random = ~ idh(Condition):animal+ Diet + ID,
                              rcov = ~idh(Condition):units,
                              pedigree = tree_drosophilidae,
                              prior = prior_rep,
                              data = data_abundances,
                              nitt = 130000*itt,
                              thin = 20*itt,
                              burnin = 30000*itt,
                              pr = TRUE)
  
  save(model_NMDS3_rep, file = here("models", "ASV_NMDS3_Repeatability.RData"))
  
} else {load(here("models", "ASV_NMDS3_Repeatability.RData"))}

NMDS3_uninjured_rep <- model_NMDS3_rep$VCV[,1] / (model_NMDS3_rep$VCV[,1] + model_NMDS3_rep$VCV[,6])
NMDS3_saline_rep <- model_NMDS3_rep$VCV[,2] / (model_NMDS3_rep$VCV[,2] + model_NMDS3_rep$VCV[,7])
NMDS3_dcv_rep <- model_NMDS3_rep$VCV[,3] / (model_NMDS3_rep$VCV[,3] + model_NMDS3_rep$VCV[,8])

# ----- 3.14 NMDS4 Heritability Model --------------------------------------------

if (!file.exists(here("models", "ASV_NMDS4_Heritability.RData"))) {
  
  model_NMDS4_her <- MCMCglmm(NMDS4 ~ Condition,
                              random = ~ idh(Condition):animal + idh(Condition):Spc.Name + Diet + ID,
                              rcov = ~idh(Condition):units,
                              pedigree = tree_drosophilidae,
                              prior = prior_her,
                              data = data_abundances,
                              nitt = 130000*itt,
                              thin = 20*itt,
                              burnin = 30000*itt,
                              pr = TRUE)
  
  save(model_NMDS4_her, file = here("models", "ASV_NMDS4_Heritability.RData"))
  
} else {load(here("models", "ASV_NMDS4_Heritability.RData"))}

NMDS4_uninjured_her <- model_NMDS4_her$VCV[,1] / (model_NMDS4_her$VCV[,1] + model_NMDS4_her$VCV[,4])
NMDS4_saline_her <- model_NMDS4_her$VCV[,2] / (model_NMDS4_her$VCV[,2] + model_NMDS4_her$VCV[,5])
NMDS4_dcv_her <- model_NMDS4_her$VCV[,3] / (model_NMDS4_her$VCV[,3] + model_NMDS4_her$VCV[,6])

# ----- 3.15 NMDS4 Repeatability Model -------------------------------------------

if (!file.exists(here("models", "ASV_NMDS4_Repeatability.RData"))) {
  
  model_NMDS4_rep <- MCMCglmm(NMDS4 ~ Condition,
                              random = ~ idh(Condition):animal+ Diet + ID,
                              rcov = ~idh(Condition):units,
                              pedigree = tree_drosophilidae,
                              prior = prior_rep,
                              data = data_abundances,
                              nitt = 130000*itt,
                              thin = 20*itt,
                              burnin = 30000*itt,
                              pr = TRUE)
  
  save(model_NMDS4_rep, file = here("models", "ASV_NMDS4_Repeatability.RData"))
  
} else {load(here("models", "ASV_NMDS4_Repeatability.RData"))}

NMDS4_uninjured_rep <- model_NMDS4_rep$VCV[,1] / (model_NMDS4_rep$VCV[,1] + model_NMDS4_rep$VCV[,6])
NMDS4_saline_rep <- model_NMDS4_rep$VCV[,2] / (model_NMDS4_rep$VCV[,2] + model_NMDS4_rep$VCV[,7])
NMDS4_dcv_rep <- model_NMDS4_rep$VCV[,3] / (model_NMDS4_rep$VCV[,3] + model_NMDS4_rep$VCV[,8])

# ----- 3.16 Summarise Models --------------------------------------------------

quant_NMDSs <- bind_rows(
  summarise_MCMC("Heritability", "NMDS 1", list(NMDS1_uninjured_her, NMDS1_saline_her, NMDS1_dcv_her)),
  summarise_MCMC("Repeatability", "NMDS 1", list(NMDS1_uninjured_rep, NMDS1_saline_rep, NMDS1_dcv_rep)),
  summarise_MCMC("Heritability", "NMDS 2", list(NMDS2_uninjured_her, NMDS2_saline_her, NMDS2_dcv_her)),
  summarise_MCMC("Repeatability", "NMDS 2", list(NMDS2_uninjured_rep, NMDS2_saline_rep, NMDS2_dcv_rep)),
  summarise_MCMC("Heritability", "NMDS 3", list(NMDS3_uninjured_her, NMDS3_saline_her, NMDS3_dcv_her)),
  summarise_MCMC("Repeatability", "NMDS 3", list(NMDS3_uninjured_rep, NMDS3_saline_rep, NMDS3_dcv_rep)),
  summarise_MCMC("Heritability", "NMDS 4", list(NMDS4_uninjured_her, NMDS4_saline_her, NMDS4_dcv_her)),
  summarise_MCMC("Repeatability", "NMDS 4", list(NMDS4_uninjured_rep, NMDS4_saline_rep, NMDS4_dcv_rep))
)

# ------------------------------------------------------------------------------
# ----- 4. Heritability/Repeatability Plots ------------------------------------
# ------------------------------------------------------------------------------

# ----- 4.1 Theme --------------------------------------------------------------

theme_herrep <- theme(aspect.ratio = 0.5,
                      text = element_text(size = 8, color = "#2e3440"),
                      axis.title = element_blank(),
                      panel.grid = element_blank(),
                      legend.title = element_text(size = 7),
                      legend.text = element_text(size = 5),
                      legend.key.size = unit(0.3, "cm"),
                      legend.position = "none",
                      strip.background = element_rect(fill = "#e5e9f0"))

# ----- 4.2 NMDS1 Plot -----------------------------------------------------------

p_model_NMDS1 <- ggplot(filter(quant_NMDSs, Metric == "NMDS 1")) +
  geom_point(aes(x = Mean, y = Condition), size = 1) +
  geom_errorbar(aes(x = Mean, xmin = HPD_Low, xmax = HPD_High, y = Condition), width = 0, size = 0.5) +
  geom_vline(aes(xintercept = 0.05), linetype = "dotted", size = 0.5) +
  facet_grid(cols = vars(Metric), rows = vars(Quant)) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme_herrep +
  guides(color = guide_legend(reverse = TRUE)) +
  ggtitle("B.")

# ----- 4.3 NMDS2 Plot -----------------------------------------------------------

p_model_NMDS2 <- ggplot(filter(quant_NMDSs, Metric == "NMDS 2")) +
  geom_point(aes(x = Mean, y = Condition), size = 1) +
  geom_errorbar(aes(x = Mean, xmin = HPD_Low, xmax = HPD_High, y = Condition), width = 0, size = 0.5) +
  geom_vline(aes(xintercept = 0.05), linetype = "dotted", size = 0.5) +
  facet_grid(cols = vars(Metric), rows = vars(Quant)) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme_herrep +
  theme(axis.text.y = element_blank()) +
  guides(color = guide_legend(reverse = TRUE)) +
  ggtitle("")

# ----- 4.4 NMDS3 Plot -----------------------------------------------------------

p_model_NMDS3 <- ggplot(filter(quant_NMDSs, Metric == "NMDS 3")) +
  geom_point(aes(x = Mean, y = Condition), size = 1) +
  geom_errorbar(aes(x = Mean, xmin = HPD_Low, xmax = HPD_High, y = Condition), width = 0, size = 0.5) +
  geom_vline(aes(xintercept = 0.05), linetype = "dotted", size = 0.5) +
  facet_grid(cols = vars(Metric), rows = vars(Quant)) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme_herrep +
  theme(axis.text.y = element_blank()) +
  guides(color = guide_legend(reverse = TRUE)) +
  ggtitle("")

# ----- 4.5 NMDS4 Plot -----------------------------------------------------------

p_model_NMDS4 <- ggplot(filter(quant_NMDSs, Metric == "NMDS 4")) +
  geom_point(aes(x = Mean, y = Condition), size = 1) +
  geom_errorbar(aes(x = Mean, xmin = HPD_Low, xmax = HPD_High, y = Condition), width = 0, size = 0.5) +
  geom_vline(aes(xintercept = 0.05), linetype = "dotted", size = 0.5) +
  facet_grid(cols = vars(Metric), rows = vars(Quant)) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme_herrep +
  theme(axis.text.y = element_blank()) +
  guides(color = guide_legend(reverse = TRUE)) +
  ggtitle("")


# ------------------------------------------------------------------------------
# ----- 5. Figure 3 Output -----------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 5.1 Arrange Figure -----------------------------------------------------

p_3A <- p_NMDS12 | p_NMDS34

p_3B <- p_model_NMDS1 | p_model_NMDS2 | p_model_NMDS3 | p_model_NMDS4

p_3B

# ----- 5.2 Save Raw Outputs ---------------------------------------------------

ggsave(here("figures", "ASV Figure 4A raw.svg"), plot = p_3A, dpi = 300, width = 7, height = 3.5)
ggsave(here("figures", "ASV Figure 4A raw.png"), plot = p_3A, dpi = 300, width = 7, height = 3.5)

ggsave(here("figures", "ASV Figure 4B raw.svg"), plot = p_3B, dpi = 300, width = 7, height = 2.5)
ggsave(here("figures", "ASV Figure 4B raw.png"), plot = p_3B, dpi = 300, width = 7, height = 2.5)


