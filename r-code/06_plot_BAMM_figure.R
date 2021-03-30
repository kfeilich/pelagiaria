# 06_plot_BAMM_figure.R

# Author: Kara Feilich
# Date modified: 20 June 2018

# Load BAMM avg rates and times 
BAMM_results <- read.csv(file="data/BAMM/rates_100.csv",header = TRUE)
labelY <- expression('Average Rate ' (lMa^{-1}))

# Build BAMM plot
bamm.plot <- ggplot(data = BAMM_results,
                          aes(x = times, y = avg, group = ID)) +
  geom_line(alpha = 0.25) +  # Nudge b/c of bin bounds
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.4)) +
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 84)) +
  scale_x_reverse(expand = c(0, 0)) +
  ggtitle("(a)")+
  labs(x='Millions of years ago', y=labelY) +
  geom_vline(xintercept = 66, colour = 'lightgrey', linetype = "dashed") +
  theme_procB() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title=element_text(face='italic'))

bamm.plot


svg(filename="figures/BAMM_AR_Depth.svg", width = 5, height = 8)
grid.draw(rbind(ggplotGrob(bamm.plot), ggplotGrob(AR_2tailed_plot), ggplotGrob(meandepth_2tailed_plot), size="last"))
dev.off()
