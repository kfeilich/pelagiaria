# Aspect ratios ###############
aspect.ratios <- read.csv ("data/PelagiariaAR_pared.csv", row.names = 1)
tree <- ape::read.nexus("data/nonconstrained_pleagia.tre")

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

tree$tip.label <- firstup(tree$tip.label) 
tree$tip.label <-  gsub('[[:digit:]]+', '', tree$tip.label)
tree$tip.label <- gsub('_$', '', tree$tip.label)
tree$tip.label <- gsub('_b$', '', tree$tip.label)

ar.names <- rownames(aspect.ratios) %>%
  .[-22]
clean.ARs <- as.data.frame(aspect.ratios[-22,])
rownames(clean.ARs) <- ar.names

aspect.ratio.formap <- treedata(phy=tree, data = clean.ARs)
aspect.ratio.named <- as.vector(aspect.ratio.formap$data)
names(aspect.ratio.named) <- rownames(aspect.ratio.formap$data)

# Get rid of duplicate tip
aspect.ratio.formap$phy <- drop.tip(aspect.ratio.formap$phy,1)

ar.map <- contMap(tree = aspect.ratio.formap$phy, x = aspect.ratio.named)
n_ar <- length(ar.map$cols)
ar.map$cols[1:n_ar] <- viridis(n_ar, direction = -1)

node_annotations <-c("Strom.","Arrip", "Bram.", "Carist.", "Chiasm.", "Icost.", "Pomat.", "Scombr.", "Scombro.", "Trichiur.", "Tetragon.", "Centro")
nodes_to_annotate <-c(133, 6, 123, 104, 95, 28, 46, 79, 64, 107, 101, 145)

pdf("figures/Aspect_ratio_contmap.pdf",width = 5,height = 8)
plot(ar.map, ftype="off",legend=FALSE,ylim=c(1-0.09*(Ntip(ar.map$tree)-1), Ntip(ar.map$tree)), xlim=c(0,1.25*max(nodeHeights(ar.map$tree))), family = "serif")
add.color.bar(40, ar.map$cols, title = "Aspect Ratio", lims = ar.map$lims, digits=1, prompt=FALSE, x=0, y = 1-0.08*(Ntip(ar.map$tree)-1), lwd=4, fsize=1, subtitle="", family="serif")
cladelabels(text = node_annotations, node = nodes_to_annotate, offset = 2,orientation = "horizontal")
dev.off()

# Import depth data #########################
depth_data <- as.data.frame(tree$tip.label) %>% 
  mutate(species = str_replace(tree$tip.label, "_", " ")) %>%
  filter(!grepl("[[:digit:]]", species)) 
depth_data <- depth_data[-1,]

# Import depth data from fishbase
depth_data_ranges <- species(depth_data$species, fields = c("Species","DepthRangeShallow", 
                                                            "DepthRangeDeep"))

depth_data <- depth_data %>%
  inner_join(depth_data_ranges, by = c("species" = "Species")) %>% 
  mutate(mean_depth = (DepthRangeShallow + DepthRangeDeep)/2,
         log_mean_depth = log10(mean_depth),
         log_max_depth = log10(DepthRangeDeep))


rownames(depth_data) <- depth_data$species
rm(depth_data_ranges)  # Clean up df after merging


# Log mean depth ###############

log.mean.depth <- select(depth_data, log_mean_depth) 
rownames(log.mean.depth) <- depth_data$"tree$tip.label"
log.mean.depth.tips <- na.exclude(log.mean.depth)

tree_fordepth <- drop.tip(tree, 1)
tree_fordepth <- extract.clade(tree_fordepth,116)
log.mean.depth.formap <- treedata(phy = tree_fordepth,data = log.mean.depth.tips)
log.mean.depth.named <- as.vector(log.mean.depth.formap$data)
names(log.mean.depth.named) <- rownames(log.mean.depth.formap$data)

mean_depth_map <- contMap(tree = log.mean.depth.formap$phy, x = log.mean.depth.named, res=200)

n_depth <- length(mean_depth_map$cols)
mean_depth_map$cols[1:n_depth] <- viridis(n_depth)[n_depth:1]

plot(mean_depth_map)

node_annotations <-c("Strom.","Arrip", "Bram.", "Carist.", "Chiasm.", "Icost.", "Pomat.", "Scombr.", "Scombro.", "Trichiur.", "Tetragon.", "Centro")
nodes_to_annotate <- c(119,5,109,90,82, 23, 36,70,56,93, 87,127)

pdf("figures/log_mean_depth_contmap.pdf",width = 5,height = 8)
plot(mean_depth_map, ftype="off",legend=FALSE,ylim=c(1-0.09*(Ntip(mean_depth_map$tree)-1), Ntip(mean_depth_map$tree)), xlim=c(0,1.25*max(nodeHeights(mean_depth_map$tree))), family = "serif")
add.color.bar(60, mean_depth_map$cols, title = "Log Mean Depth", lims = mean_depth_map$lims, digits=1, prompt=FALSE, x=0, y = 1-0.08*(Ntip(mean_depth_map$tree)-1), lwd=4, fsize=1, subtitle="", family="serif")
cladelabels(text = node_annotations, node = nodes_to_annotate, offset = 2, orientation = "horizontal")
axis(1)
dev.off()

# Log max depth ###############
log.max.depth <- select(depth_data, log_max_depth)
rownames(log.max.depth) <- depth_data$"tree$tip.label"
log.max.depth.tips <- na.exclude(log.max.depth)

log.max.depth.formap <- treedata(phy = tree, data = log.max.depth.tips)
log.max.depth.named <- as.vector(log.max.depth.formap$data)
names(log.max.depth.named) <- rownames(log.max.depth.formap$data)

contMap(tree = log.max.depth.formap$phy, x = log.max.depth.named, res = 200)

# Diet ####################
diet <- readr::read_csv("data/Pelagia_diet_data.csv") %>%
  dplyr::select(Family, Species, "Zoo-Jelly-Fish") %>%
  dplyr::rename(Zoo_Jelly_Fish = "Zoo-Jelly-Fish") %>%
  dplyr::mutate(Zoo_Jelly_Fish = as.factor(Zoo_Jelly_Fish),
                Species = str_to_title(str_extract(string = Species, pattern = "^[:alpha:]+_[:alpha:]+"))) %>%
  na.omit
diet

diet_wrownames <- as.data.frame(select(diet, Zoo_Jelly_Fish))
rownames(diet_wrownames) <- pull(diet, Species)

diet_formap <- treedata(phy = tree, data = diet_wrownames)
diet_named <- as.vector(diet_formap$data)
names(diet_named) <- rownames(diet_formap$data)

diet.simmap <- make.simmap(tree=diet_formap$phy, x = diet_named, nsim = 1000)
describe.simmap(diet.simmap, plot=TRUE)
