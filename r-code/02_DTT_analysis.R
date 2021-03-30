# 02_DTT_Analysis.R

# Author: Kara Feilich, Matt Friedman
# Date modified: 08 June 2018
# This script runs the modified MDI analysis and the rank-envelope tests of 
# disparity as defined by 
# Murrell, D. J. 2018. A global envelope test to detect non-random bursts of 
# trait evolution. Methods. Ecol. Evol. 2018;00:1–10. 
# doi: 10.1111/2041-210X.13006
# NOTE: This takes several hours to run on my laptop. 

set.seed(7)
n.simulations <- 10000  
# NB: n.simulations needs to be high or this will throw an error

# Load data ####################################################################
# load posterior sample of trees
# For unconstrained topology uncomment this ####################################
# trees <- read.nexus ("data/Posterior_100_pruned.nex")
# For constrained topology, uncomment this #####################################
trees <- read.tree ("data/pelagia_constrained_100.trees")  # For constrained topology
fix_tips <- function(tree){
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))  # For constrained topology
    x 
  }
  tree$tip.label <- firstup(tree$tip.label)
  tree$tip.label <- gsub('[[:digit:]]+', '', tree$tip.label) 
  tree$tip.label <- gsub('_$', '', tree$tip.label)  # For constrained topology
  tree$tip.label <- gsub('_b$', '', tree$tip.label)
  tree
}

trees <- lapply(trees, FUN=fix_tips)

# For all cases, continue here ######################################
trees <- lapply(trees,FUN = drop.tip, "Acanthocybium_solandri_2") 

# load aspect ratio data
aspect.ratios <- read.csv ("data/PelagiariaAR_pared.csv", row.names = 1)
rownames(aspect.ratios)[22] <- "Epinnula_affinis"  # Stand in 


# Clean all trees for Epinnula, and for matches with aspect ratio ##############
# Cleans the input trees
cleaned.trees <- clean.trees(trees, aspect.ratios)

# Run modified dtt on EACH tree in trees #######################################

# This iterates the dtt function over the list of cleaned trees, and returns the
# result as a list of each individual set of results
dtt.results.all <- purrr::map(cleaned.trees, function(x) 
  dtt1(x, aspect.ratios, plot=FALSE, nsim = n.simulations, calculateMDIp=TRUE))

# Pull the values of the MDI statistic into a vector
mdi.values <- purrr::map_dbl(dtt.results.all, getMDI)

# Pull the p-values associated with the MDI statistic into a vector
mdi.p.values <- purrr::map_dbl(dtt.results.all, getMDIPvalue)

# Bind the MDI values and the p-values into a dataframe
mdi.output <- data.frame("MDI" = mdi.values, "MDI.pval" = mdi.p.values)

# Calculate rank envelopes for each set of dtt results #########################

# This iterates the rank envelope function over the list of dtt results (the
# output from the previous function call), and returns the rank envelope result 
# as a list of each individual set of results
AR.envelope.results <- purrr::map(dtt.results.all, rank_env_dtt, 
                                    Plot = F, test = "two.sided")

# Find significant times #######################################################

# Call checkIsLower on all tree-rankenv pairs. purrr::map2 iterates the function
# (checkIsLower), over each x, y pair of (cleaned.trees, rank.envelope.results)
is.it.lower.results <- purrr::map2(cleaned.trees, AR.envelope.results,
                                   checkIsLower)
lower.runs <- purrr::map(is.it.lower.results, findLowerRuns)

is.it.higher.results <- purrr::map2(cleaned.trees, AR.envelope.results,
                                    checkIsHigher)
higher.runs <- purrr::map(is.it.higher.results, findHigherRuns)


# To get the index of the tree in cleaned.trees that produced each output, we
# need to name the elements of 'is.it.lower.results'
names(is.it.lower.results) <- as.numeric(seq(1:length(is.it.lower.results)))
names(is.it.higher.results) <- as.numeric(seq(1:length(is.it.higher.results)))

# merge all relative time points from all trees into one giant dataframe
all.lower.timepoints <- purrr::map_dfr(is.it.lower.results, rbind, .id = 'tree.index')
all.higher.timepoints <- purrr::map_dfr(is.it.higher.results, rbind, .id = 'tree.index')

# Summarize significant time points by tree for each million years #############

# First, we need maximum root age of all of the trees
unnamed.trees <- cleaned.trees  # Need to clear off the names or result is wonky
names(unnamed.trees) <- NULL
root.ages <- purrr::map_dbl(unnamed.trees, 
                             function(x) max(phytools::nodeHeights(x)))
maximum.root.age <- max(root.ages)

# Then, we need sensible bins (every million years)
all.lower.timepoints$time.bins <- cut(all.lower.timepoints$actual.time, 
                                breaks = c(0, seq(1,maximum.root.age+1)), 
                                labels=FALSE)
all.higher.timepoints$time.bins <- cut(all.higher.timepoints$actual.time, 
                                      breaks = c(0, seq(1,maximum.root.age+1)), 
                                      labels=FALSE)

# NOTE: bin labels mean that the time point was in the interval 
# (1-binlabel, bin label)

# Then, aggregate the significantly lower points by tree and by bin
agg.lower.timepts <- aggregate(x= list(is.lower = all.lower.timepoints$is.it.lower),
                         by = list(Tree = all.lower.timepoints$tree.index, 
                                   Bin = all.lower.timepoints$time.bins),
                         FUN = sum)
agg.higher.timepts <- aggregate(x= list(is.higher = all.higher.timepoints$is.it.higher),
                               by = list(Tree = all.higher.timepoints$tree.index, 
                                         Bin = all.higher.timepoints$time.bins), FUN = sum)

# Set it so multiple points from the same tree in the same bin only count as 1
agg.lower.timepts$is.lower[which(agg.lower.timepts$is.lower > 0)] <- 1   
agg.higher.timepts$is.higher[which(agg.higher.timepts$is.higher > 0)] <- 1
agg.timepts <- cbind(agg.lower.timepts, is.higher = agg.higher.timepts$is.higher)

# Get dataframe of only significant points
sig.timepts <- agg.timepts[which(agg.timepts$is.lower == 1 | agg.timepts$is.higher ==1),]
sig.timepts$Bin <- as.factor(as.character(sig.timepts$Bin))

# Aggregate by bin
sig.bins <- sig.timepts %>%
  dplyr::group_by(Bin) %>%
  dplyr::summarize(total.lower = -sum(is.lower), total.higher = sum(is.higher))%>%
  dplyr::mutate(Bin = as.integer(as.character(Bin)))

# Save out results
save(dtt.results.all, AR.envelope.results, cleaned.trees, agg.timepts, sig.bins,
     file = "AR.2tailed.rankenv_constrained.RData")

ggplot(sig.bins, aes(x=as.integer(Bin))) +
  geom_bar(aes(y=total.lower), stat= "identity", fill="black") +
  geom_bar(aes(y=total.higher), stat="identity",  fill = "black") + 
  xlab("Time Bin (million years ago)") +
  ylab("Number Posterior Trees") +
  ylim(c(-100,100))+
  scale_x_reverse(limits=c(maximum.root.age,0))+
  ggtitle("Aspect Ratio")+
  theme_custom()







