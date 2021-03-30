# 04b_Pelagiaria_w_hairtails_PCA.R

# Author:  Kara Feilich, Mimi Harrington
# Date modified: 26 November 2018
# This script conducts a principal components analysis on all taxa (including 
# hairtails), excluding those landmarks which do not apply to hairtails

lms.to.exclude <- c(4:6, 8:11)  # These LMs do not apply to hairtails

# INCLUDING THE HAIRTAILS ######################################################
# Deal with spaces in TPS IDs
pelagia.lms.wHTs <- pelagia.lms[-lms.to.exclude,,]
plotAllSpecimens(pelagia.lms.wHTs, mean = TRUE, links = NULL)  

# Find average landmarks for spp. with multiple representatives ################
# Count the number of specimens for each spp.
spp.counts.wHTs <- table(dimnames(pelagia.lms.wHTs)[[3]]) 
# Pull out the list of species for which there are more than one specimens
spp.with.mult.wHTs<- names(spp.counts.wHTs[which(spp.counts.wHTs > 1)])

# Using a for loop to average the landmarks for those spp. with mult. specimens
for (species in spp.with.mult.wHTs){  # for each spp.
  spp.average.morph.wHTs <- apply(
    pelagia.lms.wHTs[,,which(dimnames(pelagia.lms.wHTs)[[3]] == species)], 
    c(1,2), mean)  # Calculate the average landmarks
  
  # Force the averages into a 3D array with the right dimension names
  spp.average.morph.wHTs <- array(spp.average.morph.wHTs, 
                             dim = c(dim(pelagia.lms.wHTs)[1],2,1), 
                             dimnames = list(NULL, NULL, species))
  # Remove the sets of landmarks used to calculate the averages from the datset
  pelagia.lms.wHTs <- pelagia.lms.wHTs[,,- which(dimnames(pelagia.lms.wHTs)[[3]] == species)]
  # Replace the removed landmark set with the averaged one
  pelagia.lms.wHTs <- abind(pelagia.lms.wHTs, spp.average.morph.wHTs, along=3)
}
rm(species)  # Clean up after the for loop.

# Warp the dataset to account for size differences #############################
# Procrustes alignment;
pelagia.gpa.wHTs<-gpagen(pelagia.lms.wHTs, Proj = TRUE, ProcD = TRUE,
                           surfaces = NULL) 

# Visually check transformed coordinates
plotAllSpecimens(pelagia.gpa.wHTs$coords)

# Checking the deformation/Visualising the deformation
plotOutliers(pelagia.gpa.wHTs$coords)

# Perform PCA analysis #########################################################
# Output the PCA
pelagia.pca.wHTs <- plotTangentSpace(pelagia.gpa.wHTs$coords, axis1 = 1, axis2 = 2, 
                                warpgrids = TRUE, mesh = NULL, groups = NULL, 
                                verbose = FALSE , label = dimnames(pelagia.gpa.wHTs$coords)[[3]])

pelagia.pcscores.wHTs <- data.frame(pelagia.pca.wHTs$pc.scores)
pelagia.pcscores.wHTs$Taxon <- rownames(pelagia.pcscores.wHTs)
pelagia.pcscores.wHTs<- merge(pelagia.pcscores.wHTs, taxon.labels, by = "Taxon")

pelagia.pcscores.wHTs$Clade <- factor(pelagia.pcscores.wHTs$Clade,  
                                            levels = c("Scombridae", 
                                                       "Stromateoidei", 
                                                       "Trichiuroidea", 
                                                       "Other"))


# Plot ALL PCA output ##########################################################
# This code is modified from https://www.r-bloggers.com/tips-tricks-7-plotting-pca-with-tps-grids/
# Assign mean shape for use with plotRefToTarget below
mean.shape.hairtails <- mshape(pelagia.gpa.wHTs$coords)

# Make PC axis labels and limits
xlab_2 <- paste("PC 1 ", "(", round(pelagia.pca.wHTs$pc.summary$importance[2,1]*100, 1),
              "%)", sep="")
ylab_2 <- paste("PC 2 ", "(", round(pelagia.pca.wHTs$pc.summary$importance[2,2]*100, 1), 
              "%)", sep="")
xlim_2 <- c(-0.4, 0.4)
ylim_2 <- c(-0.2, 0.2)


# First item to plot: the PCA grid
plot(pelagia.pcscores.wHTs$PC1,
     pelagia.pcscores.wHTs$PC2,
     pch=21, cex=2, bg=pelagia.pcscores.wHTs$Status, xlab=xlab, ylab=ylab,
     asp=T)
legend("topleft", legend = unique(pelagia.pcscores.wHTs$Status), pch=19,  
       col=unique(pelagia.pcscores.wHTs$Status))

# Plot the warp grids

plotRefToTarget(mean.shape.hairtails, pelagia.pca.wHTs$pc.shapes$PC1min)
# Item 3
plotRefToTarget(mean.shape.hairtails,pelagia.pca.wHTs$pc.shapes$PC1max)
# Item 4
plotRefToTarget(mean.shape.hairtails, pelagia.pca.wHTs$pc.shapes$PC2min)
# Item 5
plotRefToTarget(mean.shape.hairtails, pelagia.pca.wHTs$pc.shapes$PC2max)


# Plot fossils only ############################################################
# set up shape vector
plot(pelagia.pcscores.wHTs[which(pelagia.pcscores.wHTs$Status == "fossil"),"PC1"],
     pelagia.pcscores.wHTs[which(pelagia.pcscores.wHTs$Status == "fossil"),"PC2"],
     pch=c(21,22,23,25)[factor(pelagia.pcscores.wHTs$Clade[which(
       pelagia.pcscores.wHTs$Status == "fossil")])], cex=2,
     xlab=xlab_2, ylab=ylab_2,xlim=xlim_2, ylim=ylim_2, asp=T, bg = 'red', col = 
       'white', lwd=2)
legend("topleft", 
       legend = unique(factor(pelagia.pcscores.wHTs$Clade[which(
         pelagia.pcscores.wHTs$Status == "fossil")])),
       pch = c(21,22,23,25)[unique(factor(pelagia.pcscores.wHTs$Clade[which(
         pelagia.pcscores.wHTs$Status == "fossil")]))],
       pt.bg = 'red')

# Plot modern only #############################################################
plot(pelagia.pcscores.wHTs[which(pelagia.pcscores.wHTs$Status == "modern"),"PC1"],
     pelagia.pcscores.wHTs[which(pelagia.pcscores.wHTs$Status == "modern"),"PC2"],
     pch=c(21,22,23,25)[factor(pelagia.pcscores.wHTs$Clade[which(
       pelagia.pcscores.wHTs$Status == "modern")])], cex=2,
     xlab=xlab_2, ylab=ylab_2, xlim=xlim_2, ylim=ylim_2, asp=T, bg = 'black', col = 
       'white', lwd=2)
legend("topleft", 
       legend = unique(factor(pelagia.pcscores.wHTs$Clade[which(
         pelagia.pcscores.wHTs$Status == "modern")])),
       pch = c(21,22,23,25)[unique(factor(pelagia.pcscores.wHTs$Clade[which(
         pelagia.pcscores.wHTs$Status == "modern")]))],
       pt.bg = 'black')


