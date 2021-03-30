# Model testing 

# Import trees, remove duplicate tip
trees <- read.nexus ("data/Posterior_100_pruned.nex")
trees <- lapply(trees, FUN = drop.tip, "Acanthocybium_solandri_2") 


#Import and prepare trait data for analysis
diet_colors <- c("fish" = "black", "jellyfish" = "dodgerblue", "zooplankton" = "green")

diet <- readr::read_csv("data/Pelagia_diet_data.csv") %>%
  dplyr::select(Family, Species, "Zoo-Jelly-Fish") %>%
  dplyr::rename(Zoo_Jelly_Fish = "Zoo-Jelly-Fish") %>%
  dplyr::mutate(Zoo_Jelly_Fish = as.factor(Zoo_Jelly_Fish),
                Species = str_to_title(str_extract(string = Species, pattern = "^[:alpha:]+_[:alpha:]+"))) %>%
  na.omit
diet

diet_trait <- pull(diet,Zoo_Jelly_Fish) 
names(diet_trait) <- dplyr::pull(diet, Species)  

# Match trees to trait data
for (i in c(1:length(trees))){
  trees[[i]] <- treedata(trees[[i]],data = diet_trait)$phy
  
}


# Simmap ER
er_simmaps <- map(trees, make.simmap, x = diet_trait, model = "ER", nsim = 1)

er_logL <- numeric(length = length(er_simmaps))
er_Q <- list()
for (i in c(1:length(er_simmaps))){
  er_logL[i] <- er_simmaps[[i]]$logL
  er_Q[[i]] <- er_simmaps[[i]]$Q
}

# Simmap SYM
sym_simmaps <- map(trees, make.simmap, x = diet_trait, model = "SYM", nsim = 1)

sym_logL <- numeric(length = length(sym_simmaps))
sym_Q <- list()
for (i in c(1:length(sym_simmaps))){
  sym_logL[i] <- sym_simmaps[[i]]$logL
  sym_Q[[i]] <- sym_simmaps[[i]]$Q
}


# Simmap ARD
ard_simmaps <- map(trees, make.simmap, x = diet_trait, model = "ARD", nsim = 1)
ard_logL <- numeric(length = length(ard_simmaps))
ard_Q <- list()
for (i in c(1:length(ard_simmaps))){
  ard_logL[i] <- ard_simmaps[[i]]$logL
  ard_Q[[i]] <- ard_simmaps[[i]]$Q
}

# Compile results
model_tests <- tibble(trees = trees, 
                      er_logL = er_logL, 
                      sym_logL = sym_logL, 
                      ard_logL = ard_logL,
                      er_Q = er_Q, 
                      sym_Q = sym_Q, 
                      ard_Q = ard_Q,
                      .rows = 100) %>%
  mutate(sym_er_diff = (2*(sym_logL - er_logL)) > 5.991,
         ard_er_diff = (2*(ard_logL - er_logL)) > 11.070,
         ard_sym_diff = (2*(ard_logL - sym_logL)) > 7.815)


set.seed(7)

# Runs CTT on list of trees
unnamed.trees <- trees  # Need to clear off the names or result is wonky
names(unnamed.trees) <- NULL
root.ages <- purrr::map_dbl(unnamed.trees, 
                            function(x) max(phytools::nodeHeights(x)))
maximum.root.age <- max(root.ages)
ctt_outputs <- map2(trees,c(1:length(trees)), do_ctt, diet_trait, n_sim=100)
save(ctt_outputs, file="data/ctt_outputs.RData")
