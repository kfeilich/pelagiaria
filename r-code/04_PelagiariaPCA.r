# 04_PelagiariaPCA.R

# Author: Kara Feilich, Mimi Harrington
# Date modified: 26 November 2018

# Load necessary packages #####################################################
needed.packages <- list('mvMORPH', 'geomorph', 'Morpho', 'ape', 'phytools', 
                        'surface', 'phangorn', 'ggplot2', 'ggfortify', 
                        'FactoMineR', 'Hmisc', 'data.table', 'igraph', 'dplyr', 
                        'mnormt', 'abind')
lapply(needed.packages, library, character.only = TRUE)  # loads packages 
rm(needed.packages)  # clean up

# Import taxon lookup table ####################################################
taxon.labels <- read.csv('data/Translation_table.csv', header=TRUE, 
                         colClasses = c('character', 'factor', 'factor',
                                        'factor','factor','character'))

# PCA ##########################################################################

# ANALYSIS 1: NO HAIRTAILS
# Import landmark data #########################################################
pelagia.lms <- readland.tps(file="data/__Pelagiaria_2D_landmarks_hairtail_applicable.TPS", 
                               specID=c("ID"),readcurves=FALSE, warnmsg=T)

# Remove extra spaces in TPS IDs
dimnames(pelagia.lms)[[3]] <- gsub(" ", "", 
                                              dimnames(pelagia.lms)[[3]])

# Remove the hairtails
hairtails <- taxon.labels$Taxon[which(taxon.labels$IsHairtail == 1)] # Find the hairtails
pelagia.lms.noHTs <- pelagia.lms[,,- which(dimnames(pelagia.lms)[[3]] %in% hairtails)]

# plot the raw specimens
plotAllSpecimens(pelagia.lms.noHTs, mean = TRUE, links = NULL)  

# Find average landmarks for spp. with multiple representatives ################
# Count the number of specimens for each spp.
spp.counts <- table(dimnames(pelagia.lms.noHTs)[[3]]) 
# Pull out the list of species for which there are more than one specimens
spp.with.mult <- names(spp.counts[which(spp.counts > 1)])

# Using a for loop to average the landmarks for those spp. with mult. specimens
for (species in spp.with.mult){  # for each spp.
  spp.average.morph <- apply(
    pelagia.lms.noHTs[,,which(dimnames(pelagia.lms.noHTs)[[3]] == species)], 
    c(1,2), mean)  # Calculate the average landmarks
  
  # Force the averages into a 3D array with the right dimension names
  spp.average.morph <- array(spp.average.morph, 
                             dim = c(dim(pelagia.lms.noHTs)[1],2,1), 
                             dimnames = list(NULL, NULL, species))
  # Remove the sets of landmarks used to calculate the averages from the datset
  pelagia.lms.noHTs <- pelagia.lms.noHTs[,,- which(dimnames(pelagia.lms.noHTs)[[3]] == species)]
  # Replace the removed landmark set with the averaged one
  pelagia.lms.noHTs <- abind(pelagia.lms.noHTs, spp.average.morph, along=3)
}
rm(species)  # Clean up after the for loop.

# Generalized Procrustes Analysis ##############################################
# Procrustes alignment;
pelagia.gpa.noHTs<-gpagen(pelagia.lms.noHTs, Proj = TRUE, ProcD = TRUE,
                           surfaces = NULL) 

# Visually check transformed coordinates
plotAllSpecimens(pelagia.gpa.noHTs$coords)

# Checking the deformation/Visualising the deformation
plotOutliers(pelagia.gpa.noHTs$coords)

# Perform PCA #################################################################
# Output the PCA
pelagia.pca.noHTs <- plotTangentSpace(pelagia.gpa.noHTs$coords, axis1 = 1, 
                                      axis2 = 2, warpgrids = TRUE, mesh = NULL,
                                      groups = NULL, verbose = FALSE, 
                              label = dimnames(pelagia.gpa.noHTs$coords)[[3]])

pelagia.pcscores.noHTs <- data.frame(pelagia.pca.noHTs$pc.scores)
pelagia.pcscores.noHTs$Taxon <- rownames(pelagia.pcscores.noHTs)
pelagia.pcscores.noHTs <- merge(pelagia.pcscores.noHTs, taxon.labels, by = "Taxon")
pelagia.pcscores.noHTs$Clade <- factor(pelagia.pcscores.noHTs$Clade, 
                                 levels = c("Scombridae", "Stromateoidei", 
                                            "Trichiuroidea", "Other"))

# Plot ALL PCA output ##########################################################
# This code is modified from https://www.r-bloggers.com/tips-tricks-7-plotting-pca-with-tps-grids/
# Assign mean shape for use with plotRefToTarget below
mean.shape <- mshape(pelagia.gpa.noHTs$coords)

# Make PC axis labels and limits
xlab <- paste("PC 1 ", "(", round(pelagia.pca.noHTs$pc.summary$importance[2,1]*100, 1),
              "%)", sep="")
ylab <- paste("PC 2 ", "(", round(pelagia.pca.noHTs$pc.summary$importance[2,2]*100, 1), 
              "%)", sep="")
xlim <- c(-0.3, 0.3)
ylim <- c(-0.2, 0.2)


# First item to plot: the PCA grid
plot(pelagia.pcscores.noHTs$PC1,
     pelagia.pcscores.noHTs$PC2,
     pch=21, cex=2, bg=pelagia.pcscores.noHTs$Status, xlab=xlab, ylab=ylab,
     xlim=xlim, ylim=ylim, asp=T)
legend("topleft", legend = unique(pelagia.pcscores.noHTs$Status), pch=19,  
       col=unique(pelagia.pcscores.noHTs$Status))

# Plot the warp grids

plotRefToTarget(mean.shape, pelagia.pca.noHTs$pc.shapes$PC1min)
# Item 3
plotRefToTarget(mean.shape, pelagia.pca.noHTs$pc.shapes$PC1max)
# Item 4
plotRefToTarget(mean.shape, pelagia.pca.noHTs$pc.shapes$PC2min)
# Item 5
plotRefToTarget(mean.shape, pelagia.pca.noHTs$pc.shapes$PC2max)


# Plot fossils only ############################################################
# set up shape vector
plot(pelagia.pcscores.noHTs[which(pelagia.pcscores.noHTs$Status == "fossil"),"PC1"],
     pelagia.pcscores.noHTs[which(pelagia.pcscores.noHTs$Status == "fossil"),"PC2"],
     pch=c(21,22,23,25)[factor(pelagia.pcscores.noHTs$Clade[which(
       pelagia.pcscores.noHTs$Status == "fossil")])], cex=2,
     xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, asp=T, bg = 'red', col = 
       'white', lwd=2)
legend("topleft", 
       legend = unique(factor(pelagia.pcscores.noHTs$Clade[which(
         pelagia.pcscores.noHTs$Status == "fossil")])),
  pch = c(21,22,23,25)[unique(factor(pelagia.pcscores.noHTs$Clade[which(
    pelagia.pcscores.noHTs$Status == "fossil")]))],
  pt.bg = 'red')

# Plot modern only #############################################################

plot(pelagia.pcscores.noHTs[which(pelagia.pcscores.noHTs$Status == "modern"),"PC1"],
     pelagia.pcscores.noHTs[which(pelagia.pcscores.noHTs$Status == "modern"),"PC2"],
     pch=c(21,22,23,25)[factor(pelagia.pcscores.noHTs$Clade[which(
       pelagia.pcscores.noHTs$Status == "modern")])], cex=2,
     xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, asp=T, bg = 'black', col = 
       'white', lwd=2)
legend("topleft", 
       legend = unique(factor(pelagia.pcscores.noHTs$Clade[which(
         pelagia.pcscores.noHTs$Status == "modern")])),
       pch = c(21,22,23,25)[unique(factor(pelagia.pcscores.noHTs$Clade[which(
         pelagia.pcscores.noHTs$Status == "modern")]))],
       pt.bg = 'black')




