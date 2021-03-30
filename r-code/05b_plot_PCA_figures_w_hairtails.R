#05b_plot_PCA_figures_w_hairtails

# 05_plot_PCA_figures.R

# Author: Kara Feilich
# Date modified: 12 June 2018

# Get Convex Hulls #############################################################

find_hull <- function(df) df[chull(df$PC1, df$PC2), ] 

hulls_2 <- ddply(pelagia.pcscores.wHTs, "Status", find_hull)


# Build fossils only plot ######################################################

fossil.pca.plot_2 <- ggplot(
  data = pelagia.pcscores.wHTs,
  aes(x = PC1, y = PC2, shape = Clade, fill = Clade, 
      colour=Status)) +
  geom_polygon(data=hulls_2[which(hulls_2$Status == "fossil"),], inherit.aes = FALSE,
               aes(x=PC1, y = PC2),linetype = 1, size = 0.5, fill = NA, 
               colour= "black")+
  geom_polygon(data=hulls_2[which(hulls_2$Status == "modern"),], inherit.aes = FALSE,
               aes(x=PC1, y = PC2),linetype = "dashed", size = 0.5, fill = NA, 
               colour= "grey")+
  geom_point(inherit.aes = TRUE, size=3, stroke = 0.2) +
  labs(x = xlab_2, y = ylab_2) +
  xlim(xlim_2)+
  ylim(ylim_2)+
  scale_color_manual(values = c("fossil" = "white", "modern" = NA))+
  scale_fill_manual(values = clade.colors)+
  scale_shape_manual(values = c("Scombridae" = 22, "Stromateoidei" = 23, 
                                "Trichiuroidea" = 25, "Other" = 21))+
  guides(size = FALSE, alpha = FALSE, colour = FALSE)+
  theme_procB()
fossil.pca.plot_2


# Build modern only plot ######################################################


modern.pca.plot_2 <- ggplot(
  data = pelagia.pcscores.wHTs,
  aes(x = PC1, y = PC2, shape = Clade, fill = Clade, 
      colour=Status)) +
  geom_polygon(data=hulls_2[which(hulls_2$Status == "fossil"),], inherit.aes = FALSE,
               aes(x=PC1, y = PC2),linetype = "dashed", size = 0.5, fill = NA, 
               colour= "grey")+
  geom_polygon(data=hulls_2[which(hulls_2$Status == "modern"),], inherit.aes = FALSE,
               aes(x=PC1, y = PC2),linetype = 1, size = 0.5, fill = NA, 
               colour= "black")+
  geom_point(inherit.aes = TRUE, size=3, stroke = 0.2) +
  labs(x = xlab_2, y = ylab_2) +
  xlim(xlim_2)+
  ylim(ylim_2)+
  scale_color_manual(values = c("fossil" = NA, "modern" = "white"))+
  scale_fill_manual(values = clade.colors)+
  scale_shape_manual(values = c("Scombridae" = 22, "Stromateoidei" = 23, 
                                "Trichiuroidea" = 25, "Other" = 21))+
  guides(size = FALSE, colour = FALSE, alpha=FALSE)+
  theme_procB()
modern.pca.plot_2