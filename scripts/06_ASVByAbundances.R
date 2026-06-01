# ==============================================================================
# ===== (2025) Drosophilidae 16S Virus Infection: Figure - Genera Abundances ===
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following script runs and plots individual bacterial ASV changes in
# abundance across Drosophilidae host species.

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

data_ASVs <- read_csv(here("data", "data_ASVs.csv"))

tree_drosophilidae <- read.tree(here("data", "tree_hosts.nwk"))

# ------------------------------------------------------------------------------
# ----- 1. Data Wrangling ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1 Remove non-experimental samples ------------------------------------

data_abundances <- subset_samples(phyloseq_cutoff, !(ID %in% c(as.character(307:324), "Negative")))
data_abundances <- subset_samples(data_abundances, !(Spc.Name %in% c("D. takahashii", "Z. davidi")))

# ----- 1.2. CLR Transform -----------------------------------------------------

data_presence <- as(otu_table(data_abundances), "matrix") %>% t() %>% as.data.frame()

data_presence$ID <- rownames(data_presence) %>% str_remove_all("_R1.fastq.gz") %>% as.numeric() %>% as.character()
data_presence <- gather(data_presence, key = "ASV", value = "presence", 1:(ncol(data_presence)-1))
data_presence$presence <- data_presence$presence > 0

data_clr <- as(otu_table(data_abundances), "matrix") %>% t() %>%
  cmultRepl(method = "CZM", output = "p-counts", z.warning = 1, z.delete = F)

data_clr <- clr(data_clr) %>% as.data.frame()

data_clr$ID <- rownames(data_clr) %>% str_remove_all("_R1.fastq.gz") %>% as.numeric() %>% as.character()

data_clr <- gather(data_clr, key = "ASV", value = "abundance", 1:(ncol(data_clr)-1))

data_abundances <- as.data.frame(as.matrix(sample_data(data_abundances)))

data_clr <- left_join(data_clr, data_abundances)

data_clr <- left_join(data_clr, data_presence)

data_clr$animal <- data_clr$Spc.Name %>% str_remove_all(" ")

data_clr <- left_join(data_clr, data_fly)

# ----- 1.3 Drop Tips for Species not Present ----------------------------------

unique(data_clr$animal) %in% tree_drosophilidae$tip.label

tree_drosophilidae <- drop.tip(tree_drosophilidae, setdiff(tree_drosophilidae$tip.label, data_clr$animal))

plot(tree_drosophilidae)

# ----- 1.4 Wrangle To Wide Format for MCMCglmm --------------------------------

key <- dplyr::select(data_ASVs, silva, ladderised_name)
colnames(key) <- c("ASV", "Name")
key$Name <- str_remove_all(key$Name, "\\| ")
key$Name <- str_replace_all(key$Name, " ", "_")
key$Name <- str_replace_all(key$Name, "-", "_")

data_clr_wide_presence <- select(data_clr, ID, Condition, animal, ASV, presence, Diet)
data_clr_wide_presence <- left_join(data_clr_wide_presence, key)
data_clr_wide_presence$ASV <- NULL
data_clr_wide_presence <- spread(data_clr_wide_presence, key = Name, value = presence)

data_clr_wide_abundance <- select(data_clr, ID, Condition, animal, ASV, abundance, Diet)
data_clr_wide_abundance <- left_join(data_clr_wide_abundance, key)
data_clr_wide_abundance$ASV <- NULL
data_clr_wide_abundance <- spread(data_clr_wide_abundance, key = Name, value = abundance)

data_clr <- left_join(data_clr, key)

# ------------------------------------------------------------------------------
# ----- 2. Abundance Models ----------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1 Model Prior --------------------------------------------------------

prior_pr <- list(G = list(
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000),
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000),
  G2 = list(V = diag(1), nu = 0.002)),
  R = list(V = diag(1), nu = 0.002, fix = 1))

prior_ab <- list(G = list(
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000),
  G1 = list(V = diag(1), nu = 1, alpha.mu = rep(0,1), alpha.V = diag(1) * 1000),
  G2 = list(V = diag(1), nu = 0.002)),
  R = list(V = diag(1), nu = 0.002))

# ----- 2.2 Iteration Multiplier -----------------------------------------------

itt <- 5 # Set to 10 for demonstration, 100 for publication

# ----- 2.3 Model Loop ---------------------------------------------------------

data_contrasts <- data.frame(Name = NA,
                             Condition = NA,
                             Value_pr = NA,
                             Value_pr_low = NA,
                             Value_pr_high = NA,
                             Value_ab = NA,
                             Value_ab_low = NA,
                             Value_ab_high = NA,
                             Contrast_pr = NA,
                             Contrast_pr_low = NA,
                             Contrast_pr_high = NA,
                             Contrast_pr_prob = NA,
                             Contrast_pr_prob_low = NA,
                             Contrast_pr_prob_high = NA,
                             Contrast_ab = NA,
                             Contrast_ab_low = NA,
                             Contrast_ab_high = NA,
                             Va_pr = NA,
                             Va_pr_low = NA,
                             Va_pr_high = NA,
                             Vs_pr = NA,
                             Vs_pr_low = NA,
                             Vs_pr_high = NA,
                             Vr_pr = NA,
                             Vr_pr_low = NA,
                             Vr_pr_high = NA,
                             Va_ab = NA,
                             Va_ab_low = NA,
                             Va_ab_high = NA,
                             Vs_ab = NA,
                             Vs_ab_low = NA,
                             Vs_ab_high = NA,
                             Vr_ab = NA,
                             Vr_ab_low = NA,
                             Vr_ab_high = NA)[-1,]

for(V in unique(data_clr$Name)){
  
  fixedeffs <- as.formula(paste(V, "~ Condition"))
  
  data_clr_wide_presence$Condition <- factor(data_clr_wide_presence$Condition, levels = c("None", "Ringers", "DCV"))
  data_clr_wide_abundance$Condition <- factor(data_clr_wide_abundance$Condition, levels = c("None", "Ringers", "DCV"))
  
  data_clr_wide_abundance_notzero <- data_clr_wide_abundance[data_clr_wide_presence[[V]] == TRUE, ]
  
  data_clr_wide_presence$spc <- data_clr_wide_presence$animal
  data_clr_wide_abundance_notzero$spc <- data_clr_wide_abundance_notzero$animal
  
  if (!file.exists(here("models", sprintf("%s_Contrast_Uninjected_Presence.Rdata", V)))) {
    
    model_pr <- MCMCglmm(fixed = fixedeffs,
                         random = ~animal + spc + Diet,
                         rcov = ~units,
                         pedigree = tree_drosophilidae,
                         prior = prior_pr,
                         data = data_clr_wide_presence,
                         family = "threshold",
                         nitt = 130000*itt,
                         thin = 10*itt,
                         burnin = 30000*itt,
                         pr = TRUE)
    
    save(model_pr, file = here("models", sprintf("%s_Contrast_Uninjected_Presence.Rdata", V)))
    
  } else {load(here("models", sprintf("%s_Contrast_Uninjected_Presence.Rdata", V)))}
  
  Pr_none <- as.mcmc(model_pr$Sol[,1])
  Pr_ring <- as.mcmc(model_pr$Sol[,1] + model_pr$Sol[,2])
  Pr_dcv <- as.mcmc(model_pr$Sol[,1] + model_pr$Sol[,3])
  
  Contrast_pr_ring <- as.mcmc(model_pr$Sol[,2])
  
  if (!file.exists(here("models", sprintf("%s_Contrast_Uninjected_Abundance.Rdata", V)))) {
    
    model_ab <- MCMCglmm(fixed = fixedeffs,
                         random = ~animal + spc + Diet,
                         rcov = ~units,
                         pedigree = tree_drosophilidae,
                         prior = prior_ab,
                         data = data_clr_wide_abundance_notzero,
                         family = "gaussian",
                         nitt = 130000*itt,
                         thin = 10*itt,
                         burnin = 30000*itt,
                         pr = TRUE)
    
    save(model_ab, file = here("models", sprintf("%s_Contrast_Uninjected_Abundance.Rdata", V)))
    
  } else {load(here("models", sprintf("%s_Contrast_Uninjected_Abundance.Rdata", V)))}
  
  Ab_none <- as.mcmc(model_ab$Sol[,1])
  Ab_ring <- as.mcmc(model_ab$Sol[,1] + model_ab$Sol[,2])
  Ab_dcv <- as.mcmc(model_ab$Sol[,1] + model_ab$Sol[,3])
  
  Contrast_ab_ring <- as.mcmc(model_ab$Sol[,2])
  
  data_clr_wide_presence$Condition <- factor(data_clr_wide_presence$Condition, levels = c("Ringers", "None", "DCV"))
  data_clr_wide_abundance_notzero$Condition <- factor(data_clr_wide_abundance_notzero$Condition, levels = c("Ringers", "None", "DCV"))
  
  if (!file.exists(here("models", sprintf("%s_Contrast_Saline_Presence.Rdata", V)))) {
    
    model_pr <- MCMCglmm(fixed = fixedeffs,
                         random = ~animal + spc + Diet,
                         rcov = ~units,
                         pedigree = tree_drosophilidae,
                         prior = prior_pr,
                         data = data_clr_wide_presence,
                         family = "threshold",
                         nitt = 130000*itt,
                         thin = 10*itt,
                         burnin = 30000*itt,
                         pr = TRUE)
    
    save(model_pr, file = here("models", sprintf("%s_Contrast_Saline_Presence.Rdata", V)))
    
  } else {load(here("models", sprintf("%s_Contrast_Saline_Presence.Rdata", V)))}
  
  Contrast_pr_dcv <- as.mcmc(model_pr$Sol[,3])
  
  if (!file.exists(here("models", sprintf("%s_Contrast_Saline_Abundance.Rdata", V)))) {
    
    model_ab <- MCMCglmm(fixed = fixedeffs,
                         random = ~animal + spc + Diet,
                         rcov = ~units,
                         pedigree = tree_drosophilidae,
                         prior = prior_ab,
                         data = data_clr_wide_abundance_notzero,
                         family = "gaussian",
                         nitt = 130000*itt,
                         thin = 10*itt,
                         burnin = 30000*itt,
                         pr = TRUE)
    
    save(model_ab, file = here("models", sprintf("%s_Contrast_Saline_Abundance.Rdata", V)))
    
  } else {load(here("models", sprintf("%s_Contrast_Saline_Abundance.Rdata", V)))}
  
  Contrast_ab_dcv <- as.mcmc(model_ab$Sol[,3])
  
  Contrast_pr_ring_probscale <- pnorm(Pr_ring) - pnorm(Pr_none)
  Contrast_pr_dcv_probscale <- pnorm(Pr_dcv) - pnorm(Pr_ring)
  
  V_pr <- model_pr$VCV[,1] + model_pr$VCV[,2] + model_pr$VCV[,4]
  
  Va_pr <- model_pr$VCV[,1] / V_pr
  Vs_pr <- model_pr$VCV[,2] / V_pr
  Vr_pr <- model_pr$VCV[,4] / V_pr

  V_ab <- model_ab$VCV[,1] + model_ab$VCV[,2] + model_ab$VCV[,4]
  
  Va_ab <- model_ab$VCV[,1] / V_ab
  Vs_ab <- model_ab$VCV[,2] / V_ab
  Vr_ab <- model_ab$VCV[,4] / V_ab
  
  data_contrasts <- rbind(data_contrasts,
                          data.frame(Name = V,
                                     Condition = c("None", "Ringers", "DCV"),
                                     Value_pr = c(mean(Pr_none), mean(Pr_ring), mean(Pr_dcv)),
                                     Value_pr_low = c(HPDinterval(Pr_none)[1], HPDinterval(Pr_ring)[1], HPDinterval(Pr_dcv)[1]),
                                     Value_pr_high = c(HPDinterval(Pr_none)[2], HPDinterval(Pr_ring)[2], HPDinterval(Pr_dcv)[2]),
                                     Value_ab = c(mean(Ab_none), mean(Ab_ring), mean(Ab_dcv)),
                                     Value_ab_low = c(HPDinterval(Ab_none)[1], HPDinterval(Ab_ring)[1], HPDinterval(Ab_dcv)[1]),
                                     Value_ab_high = c(HPDinterval(Ab_none)[2], HPDinterval(Ab_ring)[2], HPDinterval(Ab_dcv)[2]),
                                     Contrast_pr = c(0, mean(Contrast_pr_ring), mean(Contrast_pr_dcv)),
                                     Contrast_pr_low = c(0, HPDinterval(Contrast_pr_ring)[1], HPDinterval(Contrast_pr_dcv)[1]),
                                     Contrast_pr_high = c(0, HPDinterval(Contrast_pr_ring)[2], HPDinterval(Contrast_pr_dcv)[2]),
                                     Contrast_pr_prob = c(0, mean(Contrast_pr_ring_probscale), mean(Contrast_pr_dcv_probscale)),
                                     Contrast_pr_prob_low = c(0, HPDinterval(Contrast_pr_ring_probscale)[1], HPDinterval(Contrast_pr_dcv_probscale)[1]),
                                     Contrast_pr_prob_high = c(0, HPDinterval(Contrast_pr_ring_probscale)[2], HPDinterval(Contrast_pr_dcv_probscale)[2]),
                                     Contrast_ab = c(0, mean(Contrast_ab_ring), mean(Contrast_ab_dcv)),
                                     Contrast_ab_low = c(0, HPDinterval(Contrast_ab_ring)[1], HPDinterval(Contrast_ab_dcv)[1]),
                                     Contrast_ab_high = c(0, HPDinterval(Contrast_ab_ring)[2], HPDinterval(Contrast_ab_dcv)[2]),
                                     Va_pr = mean(Va_pr),
                                     Va_pr_low = HPDinterval(Va_pr)[1],
                                     Va_pr_high = HPDinterval(Va_pr)[2],
                                     Vs_pr = mean(Vs_pr),
                                     Vs_pr_low = HPDinterval(Vs_pr)[1],
                                     Vs_pr_high = HPDinterval(Vs_pr)[2],
                                     Vr_pr = mean(Vr_pr),
                                     Vr_pr_low = HPDinterval(Vr_pr)[1],
                                     Vr_pr_high = HPDinterval(Vr_pr)[2],
                                     Va_ab = mean(Va_ab),
                                     Va_ab_low = HPDinterval(Va_ab)[1],
                                     Va_ab_high = HPDinterval(Va_ab)[2],
                                     Vs_ab = mean(Vs_ab),
                                     Vs_ab_low = HPDinterval(Vs_ab)[1],
                                     Vs_ab_high = HPDinterval(Vs_ab)[2],
                                     Vr_ab = mean(Vr_ab),
                                     Vr_ab_low = HPDinterval(Vr_ab)[1],
                                     Vr_ab_high = HPDinterval(Vr_ab)[2]))
  
}

# ------------------------------------------------------------------------------
# ----- 3. Figure 5-6 ----------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 3.1. Wrangle Plot Elements ---------------------------------------------

data_contrasts$Condition <- factor(data_contrasts$Condition,
                                   levels = c("None", "Ringers", "DCV"))

data_contrasts$Name <- factor(data_contrasts$Name, levels = rev(sort(unique(data_contrasts$Name))))

# ----- 3.2. Detection Probability Plot ----------------------------------------

p1 <- ggplot(filter(data_contrasts, Condition == "None")) +
  geom_errorbar(aes(xmin = pnorm(Value_pr_low), xmax = pnorm(Value_pr_high), y = Name), width = 0, color = "#2c3e50") +
  geom_point(aes(x = pnorm(Value_pr), y = Name), color = "#2c3e50", size = 1) +
  scale_x_continuous(name = "Detection Probability") +
  facet_wrap(~Condition) +
  theme_bw() +
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 8, color = "#2e3440"),
        strip.background = element_rect(fill = "#e5e9f0"))

# ----- 3.3. Change in Detection Probability Plots -----------------------------

p2 <- ggplot(filter(data_contrasts, Condition == "Ringers")) +
  geom_errorbar(aes(xmin = Contrast_pr_prob_low, xmax = Contrast_pr_prob_high, y = Name), width = 0, color = "#bdc3c7") +
  geom_point(aes(x = Contrast_pr_prob, y = Name), color = "#bdc3c7", size = 1) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(name = "Change in Detection Probability", limits = c(-0.4, 0.4), expand = c(0,0), breaks = c(-0.25, 0, 0.25)) +
  facet_wrap(~Condition) +
  theme_bw() +
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 8, color = "#2e3440"),
        strip.background = element_rect(fill = "#e5e9f0"))

p3 <- ggplot(filter(data_contrasts, Condition == "DCV")) +
  geom_errorbar(aes(xmin = Contrast_pr_prob_low, xmax = Contrast_pr_prob_high, y = Name), width = 0, color = "#bdc3c7") +
  geom_point(aes(x = Contrast_pr_prob, y = Name), color = "#bdc3c7", size = 1) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(name = "Change in Detection Probability", limits = c(-0.4, 0.4), expand = c(0,0), breaks = c(-0.25, 0, 0.25)) +
  scale_y_discrete(position = "right") +
  facet_wrap(~Condition) +
  theme_bw() +
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"),
        text = element_text(size = 8, color = "#2e3440"),
        strip.background = element_rect(fill = "#e5e9f0"))

# ----- 3.4. Abundance Plot ----------------------------------------------------

p4 <- ggplot(filter(data_contrasts, Condition == "None")) +
  geom_errorbar(aes(xmin = Value_ab_low, xmax = Value_ab_high, y = Name), width = 0, color = "#2c3e50") +
  geom_point(aes(x = Value_ab, y = Name), color = "#2c3e50", size = 1) +
  scale_x_continuous(name = "CLR-transformed Abundance") +
  facet_wrap(~Condition) +
  theme_bw() +
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 8, color = "#2e3440"),
        strip.background = element_rect(fill = "#e5e9f0"))

# ----- 3.4. Change in Abundance Plots -----------------------------------------

p5 <- ggplot(filter(data_contrasts, Condition == "Ringers")) +
  geom_errorbar(aes(xmin = Contrast_ab_low, xmax = Contrast_ab_high, y = Name), width = 0, color = "#bdc3c7") +
  geom_point(aes(x = Contrast_ab, y = Name), color = "#bdc3c7", size = 1) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(name = "Change in CLR-transformed Abundance", limits = c(-4.25, 4.25), expand = c(0,0)) +
  facet_wrap(~Condition) +
  theme_bw() +
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 8, color = "#2e3440"),
        strip.background = element_rect(fill = "#e5e9f0"))

p6 <- ggplot(filter(data_contrasts, Condition == "DCV")) +
  geom_errorbar(aes(xmin = Contrast_ab_low, xmax = Contrast_ab_high, y = Name), width = 0, color = "#bdc3c7") +
  geom_point(aes(x = Contrast_ab, y = Name), color = "#bdc3c7", size = 1) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_x_continuous(name = "Change in CLR-transformed Abundance", limits = c(-5, 5), expand = c(0,0)) +
  scale_y_discrete(position = "right") +
  facet_wrap(~Condition) +
  coord_cartesian(xlim = c(-4.25, 4.25)) +
  theme_bw() +
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"),
        text = element_text(size = 8, color = "#2e3440"),
        strip.background = element_rect(fill = "#e5e9f0"))


# ----- 3.5 Save Raw Outputs ---------------------------------------------------

plot <- (p1 | p2 | p3) / (p4 | p5 | p6)

plot

ggsave(here("figures", "ASV Figure 5_6 raw.svg"), plot = plot, dpi = 300, width = 7, height = 9)
ggsave(here("figures", "ASV Figure 5_6 raw.png"), plot = plot, dpi = 300, width = 7, height = 9)