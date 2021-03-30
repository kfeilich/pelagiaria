# 06_plot_composite_figure.R

# Author: Kara Feilich
# Date modified: 12 June 2018

blank.panel <- rectGrob(gp=gpar(col=NA))

# Get warp grids
ref <- mshape(pelagia.gpa.noHTs$coords)

svg("figures/far_left_Aanzac.svg") 
plotRefToTarget(ref, pelagia.gpa.noHTs$coords[,,"Assurger_anzac"])
dev.off()

svg("figures/top_left_Btenuis.svg")
plotRefToTarget(ref, pelagia.gpa.noHTs$coords[,,"Benthodesmus_tenuis"])
dev.off()

svg("figures/top_right_Ppetersii.svg")
plotRefToTarget(ref, pelagia.gpa.noHTs$coords[,,"Pterycombus_petersii"])
dev.off()

svg("figures/bottom_right_Bcaribbea.svg")
plotRefToTarget(ref, pelagia.gpa.noHTs$coords[,,"Brama_caribbea"])
dev.off()

svg("figures/bottom_Pelamensis.svg")
plotRefToTarget(ref, pelagia.gpa.noHTs$coords[,,"Paucaichthys_elamensis"])
dev.off()

svg("figures/bottom_left_Ecancellata.svg")
plotRefToTarget(ref, pelagia.gpa.noHTs$coords[,,"Epinnula_cancellata"])
dev.off()

svg("figures/bottom_modern_Tatlanticus.svg")
plotRefToTarget(ref, pelagia.gpa.noHTs$coords[,,"Tetragonurus_atlanticus"])
dev.off()


# remove legends
fossil.pca.panel <- fossil.pca.plot + 
  guides(shape=FALSE,color=FALSE, fill = FALSE) +
  annotate(geom = 'text', x = 0.18, y = -0.17, label = 'Paleocene-Eocene', 
           family = 'TimesNewRoman')
#  theme(legend.justification = c(0,0), legend.position = c(0.01,0.01))
modern.pca.panel <- modern.pca.plot +
  guides(shape=FALSE,color=FALSE, fill = FALSE)+
  annotate(geom = 'text', x = 0.2, y = -0.17, label = 'Recent', 
           family = 'TimesNewRoman') # + 
  # annotation_custom(far_right, xmin=0.3)
  # annotation_custom(bottomright, xmin=0.1, xmax = 0.45,ymin =-0.3, ymax=-0.05)+
  # annotation_custom(topleft, xmin=-0.55, xmax = -0.2, ymin=-0.02, ymax=0.23)+
  # annotation_custom(topright, xmin=0.1, xmax = 0.45,ymin =-0.02, ymax=0.23)
  
timepoints.panel <- timepoints.plot + 
  annotate(geom = 'text', x = 61, y = 0.75, label = 'K-Pg', 
           family = 'TimesNewRoman')
  

# build composite
svg("figures/Composite_BAMM_morphospace_nohairtails.svg", width = 7, height=4)
grid.arrange(bamm.plot, timepoints.panel, fossil.pca.panel, modern.pca.panel,
             nrow=2, heights = c(2/5, 3/5)) 
dev.off()


# Same, for hairtails
ref_2 <- mshape(pelagia.gpa.wHTs$coords)
PC_2 <- pelagia.pca.wHTs$pc.scores[,1:2]
preds_2 <- shape.predictor(pelagia.gpa.wHTs$coords, x = PC_2, Intercept = FALSE,
                         bottomleft = c(min(PC_2[,1]), min(PC_2[,2])),
                         bottomright = c(max(PC_2[,1]), min(PC_2[,2])),
                         topleft = c(min(PC_2[,1]), max(PC_2[,2])),
                         topright = c(max(PC_2[,1]), max(PC_2[,2])))
bottomleft_2 <- base2grob(function() plotRefToTarget(ref_2, preds_2$bottomleft))
bottomright_2 <- base2grob(function() plotRefToTarget(ref_2, preds_2$bottomright))
topleft_2 <- base2grob(function() plotRefToTarget(ref_2, preds_2$topleft))
topright_2 <- base2grob(function() plotRefToTarget(ref_2, preds_2$topright))

fossil.pca.panel_2 <- fossil.pca.plot_2 +  guides(shape=FALSE,color=FALSE, fill = FALSE) # theme(legend.justification = c(0,0), legend.position = c(0.01,0.01))
modern.pca.panel_2 <- modern.pca.plot_2 + guides(shape=FALSE,color=FALSE, fill=FALSE)  # + 
  # annotation_custom(bottomleft_2, xmin=-0.55, xmax = -0.2, ymin=-0.3, ymax=-0.05)+
  # annotation_custom(bottomright_2, xmin=0.1, xmax = 0.45,ymin =-0.3, ymax=-0.05)+
  # nnotation_custom(topleft_2, xmin=-0.55, xmax = -0.2, ymin=-0.02, ymax=0.23)+
  # annotation_custom(topright_2, xmin=0.1, xmax = 0.45,ymin =-0.02, ymax=0.23)

svg("figures/Morphospaces_hairtails.svg", width = 10, height= 4)
grid.arrange(fossil.pca.panel_2, modern.pca.panel_2,
             nrow=1) 
dev.off()

fossil.pca.panel_3 <- fossil.pca.plot_2
svg("figures/Morphospaces_hairtails_forlegend.svg", width = 10, height= 4)
grid.arrange(fossil.pca.panel_3, modern.pca.panel_2,
             nrow=1) 
dev.off()


ref2 <- mshape(pelagia.gpa.wHTs$coords)

svg("figures/HT_farright_Pargenteus.svg")
plotRefToTarget(ref2, pelagia.gpa.wHTs$coords[,,"Pampus_argenteus"])
dev.off()

svg("figures/HT_farleft_Aanzac.svg")
plotRefToTarget(ref2, pelagia.gpa.wHTs$coords[,,"Assurger_anzac"])
dev.off()

svg("figures/HT_topmid_Cmacropus.svg")
plotRefToTarget(ref2, pelagia.gpa.wHTs$coords[,,"Caristius_macropus"])
dev.off()

svg("figures/HT_bottom_Ecancellata.svg")
plotRefToTarget(ref2, pelagia.gpa.wHTs$coords[,,"Epinnula_cancellata"])
dev.off()

