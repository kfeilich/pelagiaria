# 02a_depth_analysis.R

# Author: Kara Feilich
# Date modified: 19 February 2019
# This script runs the modified MDI analysis and the rank-envelope tests of
# depth reconstructions as defined by 
# Murrell, D. J. 2018. A global envelope test to detect non-random bursts of 
# trait evolution. Methods. Ecol. Evol. 2018;00:1–10. 
# doi: 10.1111/2041-210X.13006
# NOTE: This takes several hours to run on my laptop. 

set.seed(7)
n.simulations <- 10000
# NB: n.simulations needs to be high or this will throw an error

# Load data ####################################################################
# load posterior sample of trees
trees <- read.nexus ("data/Posterior_100_pruned.nex") 
trees <- lapply(trees,FUN = drop.tip, "Acanthocybium_solandri_2")  # Remove duplicate species
tips <- trees[[1]]$tip.label

# Set up data frame, excluding species with numbers.
depth_data <- as.data.frame(tips) %>% 
  mutate(species = str_replace(tips, "_", " ")) %>%
  filter(!grepl("[[:digit:]]", species)) 


# Import depth data from fishbase
depth_data_ranges <- species(depth_data$species, fields = c("Species", "DepthRangeShallow", 
                                                            "DepthRangeDeep"))

depth_data <- depth_data %>%
  inner_join(depth_data_ranges, by = c("species" = "Species")) %>% 
  mutate(mean_depth = (DepthRangeShallow + DepthRangeDeep)/2,
         log_mean_depth = log10(mean_depth),
         log_max_depth = log10(DepthRangeDeep))

rm(depth_data_ranges)  # Clean up df after merging

# Plot ditribution of mean depths through time. 

# Make single column dataframes with names
mean.depth <- dplyr::select(depth_data, mean_depth)
rownames(mean.depth) <- depth_data$tip.label

log.mean.depth <- dplyr::select(depth_data, log_mean_depth)
rownames(log.mean.depth) <- depth_data$tip.label

max.depth <- dplyr::select(depth_data, DepthRangeDeep)
rownames(max.depth) <- depth_data$tip.label

log.max.depth <- dplyr::select(depth_data, log_max_depth)
rownames(log.max.depth) <- depth_data$tip.label

mean_depth_dtt_twoside <- run_modified_dtt(mean.depth, trees, 
                                           n_sim = n.simulations, 
                                           test_type = "two.sided")
save(mean_depth_dtt_twoside, file = "mean_depth_twoside_res.RData")
rm(mean_depth_dtt_twoside)

log_mean_depth_dtt_twoside <- run_modified_dtt(log.mean.depth, trees, 
                                               n_sim = n.simulations, 
                                               test_type = "two.sided")
save(log_mean_depth_dtt_twoside, file = "log_mean_depth_twoside_res.RData")
rm(log_mean_depth_dtt_twoside)

max_depth_dtt_twoside <- run_modified_dtt(max.depth, trees, 
                                          n_sim = n.simulations, 
                                          test_type = "two.sided")
save(max_depth_dtt_twoside, file = "max_depth_twoside_res.RData")
rm(max_depth_dtt_twoside)

log_max_depth_dtt_twoside <- run_modified_dtt(log.max.depth, trees, 
                                              n_sim = n.simulations, 
                                              test_type = "two.sided")
save(log_max_depth_dtt_twoside, file = "log_max_depth_twoside_res.RData")
rm(log_max_depth_dtt_twoside)


# Do without apparent outlier taxa

trees.noShallow <- vector("list", length=100)
for (i in c(1:length(trees))) {
  trees.noShallow[[i]] <- drop.tip(trees[[i]],
                                   tip = c("Scomberomorus_maculatus", 
                                           "Scomberomorus_regalis", 
                                           "Acanthocybium_solandri",
                                           "Paradiplospinus_antarcticus", 
                                           "Aphanops_carbo"))
}

class(trees.noShallow) <- "multiPhylo"

log_mean_depth_dtt_twoside_noshallow <- run_modified_dtt(log.mean.depth, 
                                                         trees.noShallow, 
                                                         n_sim = n.simulations, 
                                                         test_type = "two.sided")
save(log_mean_depth_dtt_twoside_noshallow,
     file = "log_mean_depth_twoside_res_noshallow.RData")
