# 05_plot_PCA_figures.R

# Author: Kara Feilich
# Date modified: 12 June 2018

# Get Convex Hulls #############################################################

find_hull <- function(df) df[chull(df$PC1, df$PC2), ] 

hulls <- ddply(pelagia.pcscores.noHTs, "Status", find_hull)
hulls2 <- ddply(pelagia.pcscores.noHTs, "InTimeWindow", find_hull)

# Make color palette ##############################
clade.colors <- c("Scombridae" = '#E69F00', "Stromateoidei" = "#56B4E9", "Trichiuroidea" = "#009E73", "Other" = "#F0E442")

# Build fossils only plot ######################################################

fossil.pca.plot <- ggplot(
  data = pelagia.pcscores.noHTs,
  aes(x = PC1, y = PC2, shape = Clade, fill = Clade, 
      colour=Status)) +
  geom_polygon(data=hulls[which(hulls$Status == "fossil"),], inherit.aes = FALSE,
               aes(x=PC1, y = PC2),linetype = 1, size = 0.5, fill = NA, 
               colour= "black")+
  geom_polygon(data=hulls[which(hulls$Status == "modern"),], inherit.aes = FALSE,
               aes(x=PC1, y = PC2),linetype = 'dashed', size = 0.5, fill = NA, 
               colour= "grey")+
  geom_point(inherit.aes = TRUE, size = 3, stroke = 0.2) +
  labs(x = xlab, y = ylab) +
  xlim(xlim)+
  ylim(ylim)+
  scale_color_manual(values = c("fossil" = "white", "modern" = NA))+
  scale_fill_manual(values = clade.colors)+
  scale_shape_manual(values = c("Scombridae" = 22, "Stromateoidei" = 23, 
                                "Trichiuroidea" = 25, "Other" = 21))+
  guides(size = FALSE, alpha=FALSE,  colour=FALSE)+
  theme_procB()
fossil.pca.plot 

# Build modern only plot ######################################################


modern.pca.plot <- ggplot(
  data = pelagia.pcscores.noHTs,
  aes(x = PC1, y = PC2, shape = Clade, fill = Clade, 
      colour=Status)) +
  geom_polygon(data=hulls[which(hulls$Status == "fossil"),], inherit.aes = FALSE,
               aes(x=PC1, y = PC2),linetype = 'dashed', size = 0.5, fill = NA, 
               colour= "grey")+
  geom_polygon(data=hulls[which(hulls$Status == "modern"),], inherit.aes = FALSE,
               aes(x=PC1, y = PC2),linetype = 1, size = 0.5, fill = NA, 
               colour= "black")+
  geom_point(inherit.aes = TRUE, size = 3, stroke = 0.2) +
  labs(x = xlab, y = ylab) +
  xlim(xlim)+
  ylim(ylim)+
  scale_color_manual(values = c("fossil" = NA, "modern" = "white"))+
  scale_fill_manual(values = clade.colors)+
  scale_shape_manual(values = c("Scombridae" = 22, "Stromateoidei" = 23, 
                                "Trichiuroidea" = 25, "Other" = 21))+
  guides(size = FALSE, alpha=FALSE, color = FALSE)+
  theme_procB()

modern.pca.plot
# Plot all specimens ###########################################################
all.pca.plot <- ggplot(
  data = pelagia.pcscores.noHTs,
  aes(x = PC1, y = PC2, shape = Clade, fill = Clade, 
      colour=Status)) +
  geom_polygon(data=hulls[which(hulls$Status == "fossil"),], inherit.aes = FALSE,
               aes(x=PC1, y = PC2),linetype = "dashed", size = 0.5, fill = NA, 
               colour= "grey")+
  geom_polygon(data=hulls[which(hulls$Status == "modern"),], inherit.aes = FALSE,
               aes(x=PC1, y = PC2),linetype = 1, size = 0.5, fill = NA, 
               colour= "black")+
  geom_point(inherit.aes = TRUE, size = 3, stroke = 0.4) +
  labs(x = xlab, y = ylab) +
  xlim(xlim)+
  ylim(ylim)+
  scale_color_manual(values = c("fossil" = "white", "modern" = "black"))+
  scale_fill_manual(values = clade.colors)+
  scale_shape_manual(values = c("Scombridae" = 22, "Stromateoidei" = 23, 
                                "Trichiuroidea" = 25, "Other" = 21))+
  guides(size = FALSE, alpha=FALSE, color = FALSE)+
  theme_procB()
all.pca.plot

################################################################################
timewindow.pca.plot <- ggplot(
  data = pelagia.pcscores.noHTs,
  aes(x = PC1, y = PC2, shape = Clade, fill = Clade, 
      colour=Status)) +
  geom_polygon(data=hulls2[which(hulls2$InTimeWindow == 1),], inherit.aes = FALSE,
               aes(x=PC1, y = PC2),linetype = "dashed", size = 0.5, fill = NA, 
               colour= "grey")+
  geom_polygon(data=hulls2[which(hulls2$InTimeWindow ==0),], inherit.aes = FALSE,
               aes(x=PC1, y = PC2),linetype = 1, size = 0.5, fill = NA, 
               colour= "black")+
  geom_point(inherit.aes = TRUE, size = 3, stroke = 0.4) +
  labs(x = xlab, y = ylab) +
  xlim(xlim)+
  ylim(ylim)+
  scale_color_manual(values = c("fossil" = "white", "modern" = "black"))+
  scale_fill_manual(values = clade.colors)+
  scale_shape_manual(values = c("Scombridae" = 22, "Stromateoidei" = 23, 
                                "Trichiuroidea" = 25, "Other" = 21))+
  guides(size = FALSE, alpha=FALSE, color = FALSE)+
  theme_procB()
timewindow.pca.plot