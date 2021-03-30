# Log mean depth check ##############################

load("data/log_mean_depth_twoside_res.RData")

unnamed.trees <- log_mean_depth_dtt_twoside$cleaned_trees  # Need to clear off the names or result is wonky
names(unnamed.trees) <- NULL
root.ages <- purrr::map_dbl(unnamed.trees, 
                            function(x) max(phytools::nodeHeights(x)))
maximum.root.age <- max(root.ages)

agg.timepts <- cbind(log_mean_depth_dtt_twoside$agg_lower_timepoints, 
                     is.higher = log_mean_depth_dtt_twoside$agg_higher_timepoints$is.higher)

# Get dataframe of only significant points
sig.timepts <- agg.timepts[which(agg.timepts$is.lower == 1 | agg.timepts$is.higher ==1),]
sig.bins <- sig.timepts %>%
  dplyr::group_by(Bin) %>%
  dplyr::summarize(total.lower = -sum(is.lower), total.higher = sum(is.higher)) %>%
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
zeros_lower <- data.frame(Bin = c(1:79)[which(!(c(1:79) %in% cumulative_lower$Bin))],
                          Freq = 0)
zeros_higher <- data.frame(Bin = c(1:79)[which(!(c(1:79) %in% cumulative_higher$Bin))],
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


meandepth_2tailed_plot <- ggplot(sig.bins, aes(x=Bin)) +
  geom_area(aes(x=Bin, y=cumFreq), data = cumulative_higher, fill="lightgrey")+
 geom_area(aes(x=Bin, y=cumFreq), data = cumulative_lower, fill = "lightgrey")+
   geom_bar(aes(y=-total.lower), stat= "identity", fill="black") +
  geom_bar(aes(y=total.higher), stat="identity",  fill = "black") + 
  xlab("Time Bin (million years ago)") +
  ylab("Number Posterior Trees") +
  ylim(c(-100,100))+
  scale_x_reverse(limits=c(maximum.root.age,0), expand=c(0,0))+
  geom_hline(yintercept = 0, lwd = 0.2)+
  geom_vline(xintercept = 66, colour = 'lightgrey', linetype = "dashed") +
  ggtitle("(c)")+
  theme_procB()+ 
  theme(plot.title=element_text(face='italic'))

gggeo_scale(meandepth_2tailed_plot, pos="bottom")

rm(log_mean_depth_dtt_twoside)
# Log max depth check #######################################

load("data/log_max_depth_twoside_res.RData")

agg.timepts <- cbind(log_max_depth_dtt_twoside$agg_lower_timepoints, 
                     is.higher = log_max_depth_dtt_twoside$agg_higher_timepoints$is.higher)

# Get dataframe of only significant points
# Get dataframe of only significant points
sig.timepts <- agg.timepts[which(agg.timepts$is.lower == 1 | agg.timepts$is.higher ==1),]
sig.bins <- sig.timepts %>%
  dplyr::group_by(Bin) %>%
  dplyr::summarize(total.lower = -sum(is.lower), total.higher = sum(is.higher)) %>%
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
zeros_lower <- data.frame(Bin = c(1:79)[which(!(c(1:79) %in% cumulative_lower$Bin))],
                          Freq = 0)
zeros_higher <- data.frame(Bin = c(1:79)[which(!(c(1:79) %in% cumulative_higher$Bin))],
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


maxdepth_2tailed_plot <- ggplot(sig.bins, aes(x=Bin)) +
  geom_area(aes(x=Bin, y=cumFreq), data = cumulative_higher, fill="lightgrey")+
  geom_area(aes(x=Bin, y=cumFreq), data = cumulative_lower, fill = "lightgrey")+
  geom_bar(aes(y=-total.lower), stat= "identity", fill="black") +
  geom_bar(aes(y=total.higher), stat="identity",  fill = "black") + 
  xlab("Time Bin (million years ago)") +
  ylab("Number Posterior Trees") +
  ylim(c(-100,100))+
  scale_x_reverse(limits=c(maximum.root.age,0))+
  geom_hline(yintercept = 0, lwd = 0.2)+
  geom_vline(xintercept = 66, colour = 'lightgrey', linetype = "dashed") +
  ggtitle("Log Max Depth")+
  theme_custom()

maxdepth_2tailed_plot 

rm(log_max_depth_dtt_twoside)


# Without shallow taxa

load("data/log_mean_depth_twoside_res_noshallow.RData.RData")
