# 02b_load_DTT_analysis.R

# Author: Kara Feilich
# Date modified: 12 June 2018

# load(file = "data/AR.2tailed.rankenv.RData")  # unconstrained
load(file = "AR.2tailed.rankenv_constrained.RData")  #constrained

unnamed.trees <- cleaned.trees  # Need to clear off the names or result is wonky
names(unnamed.trees) <- NULL
root.ages <- purrr::map_dbl(unnamed.trees, 
                            function(x) max(phytools::nodeHeights(x)))
maximum.root.age <- max(root.ages)

sig.bins <- sig.bins %>%
  mutate(Bin = as.integer(as.character(Bin)))


# Find the number of trees lower and higher by each bin #################
cumulative <- data.frame(Tree = c(1:100), max.lower = NA, max.higher = NA)
for (i in 1:length(cumulative$Tree)){
  cumulative$max.lower[i] <- max(agg.timepts$Bin[which(agg.timepts$Tree == i & agg.timepts$is.lower==1)])
  cumulative$max.higher[i] <- max(agg.timepts$Bin[which(agg.timepts$Tree == i & agg.timepts$is.higher==1)])
}

cumulative$max.lower[is.infinite(cumulative$max.lower)] <- NA
cumulative$max.higher[is.infinite(cumulative$max.higher)] <- NA

# Sum up the cumulative number lower or higher by each bin 

cumulative_lower <- data.frame(table(cumulative$max.lower)) ## Num trees lower than envelope initially in that bin
names(cumulative_lower) <- c("Bin", "Freq")
cumulative_lower$Bin <- as.integer(as.character(cumulative_lower$Bin))
cumulative_higher <- data.frame(table(cumulative$max.higher))  ## Num trees higher than envelope initially in that bin
names(cumulative_higher) <- c("Bin", "Freq")
cumulative_higher$Bin <- as.integer(as.character(cumulative_higher$Bin))

# Add the bins with 0s
zeros_lower <- data.frame(Bin = c(1:max(agg.timepts$Bin))[which(!(c(1:max(agg.timepts$Bin)) %in% cumulative_lower$Bin))],
                          Freq = 0)
zeros_higher <- data.frame(Bin = c(1:max(agg.timepts$Bin))[which(!(c(1:max(agg.timepts$Bin)) %in% cumulative_higher$Bin))],
                          Freq = 0)
cumulative_lower <- arrange(rbind(cumulative_lower, zeros_lower), desc(Bin))
cumulative_higher <- arrange(rbind(cumulative_higher, zeros_higher), desc(Bin))

# Find running tally
cumulative_lower$cumFreq <- NA
cumulative_higher$cumFreq <- NA
for (i in 1:nrow(cumulative_lower)){
  if(i==1){
    cumulative_lower$cumFreq[i] <-  cumulative_lower$Freq[i]
    cumulative_higher$cumFreq[i] <- cumulative_higher$Freq[i]
  } else {
    cumulative_lower$cumFreq[i] <- cumulative_lower$cumFreq[i-1] + cumulative_lower$Freq[i]
    cumulative_higher$cumFreq[i] <- cumulative_higher$cumFreq[i-1] + cumulative_higher$Freq[i]
  }
}

plot_AR_trees_tally <- ggplot()+
  geom_area(aes(x=Bin, y=cumFreq), data = cumulative_higher)+
  geom_area(aes(x=Bin, y=-cumFreq), data=cumulative_lower, fill = "lightgrey")+
  ylim(c(-100,100)) + 
  ggtitle("Aspect Ratio")+
  scale_x_reverse()+
  ylab("Cumulative Number Post. Trees")+
  theme_custom()


AR_2tailed_plot <- ggplot(sig.bins, aes(x=Bin)) +
  geom_bar(aes(y=total.lower), stat= "identity", fill="black") +
  geom_bar(aes(y=total.higher), stat="identity",  fill = "black") + 
  geom_area(aes(x=Bin, y=-cumFreq), data = cumulative_lower,alpha=0.2)+
  geom_area(aes(x=Bin, y=cumFreq), data = cumulative_higher,alpha=0.2)+
  ylab("Number Posterior Trees") +
  ylim(c(-100,100))+
  scale_x_reverse(limits=c(maximum.root.age,0), expand= c(0,0))+
  ggtitle("(b)")+
  geom_vline(xintercept = 66, colour = 'lightgrey', linetype = "dashed") +
  geom_hline(yintercept = 0, lwd = 0.2)+
  theme_procB()+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), plot.title=element_text(face='italic'))

rm(dtt.results.all, AR.envelope.results, cleaned.trees)
