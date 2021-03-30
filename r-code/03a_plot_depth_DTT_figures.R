# 03_plot_DTT_figures.R

# Author: Kara Feilich
# Date modified: 08 June 2018
# Run this after 02_DTT_analysis.R to plot interesting things

# Derive MDI values ##########################################################
# Pull the values of the MDI statistic into a vector
mdi.depth.values <- purrr::map_dbl(dtt.depth.results.all, getMDI)
# Pull the p-values associated with the MDI statistic into a vector
mdi.depth.p.values <- purrr::map_dbl(dtt.depth.results.all, getMDIPvalue)
# Bind the MDI values and the p-values into a dataframe
mdi.depth.output <- data.frame("MDI" = mdi.depth.values, "MDI.pval" = mdi.depth.p.values)

# Plot all of the significantly lower points from rank envelope#################
timepoints.depth.plot <- ggplot(data = lower.depth.bins,
                          aes(x = Bin, y = NumTrees/length(cleaned.depth.trees))) +
  geom_col(position=position_nudge(x=0.5), fill = "grey50") +  # Nudge b/c of bin bounds
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.0)) +  
  labs(x='Millions of years ago (binned)', y='Proportion of\nPosterior Trees') +
  scale_x_reverse()+
  geom_vline(xintercept = 66, colour = 'red', linetype = "dashed") +
  theme_procB()

timepoints.depth.plot

# Make rank envelope p-value range histograms ##################################

# Pulls all p-values from list of rank-envelope results into a vector (type double)
depth.p.values <- purrr::map_dbl(depth.rank.envelope.results, getPvalue)

# Pulls all p-value intervals from list of rank envelope results into a dataframe
# with columns of Lower (bound) and Upper (bound)
depth.p.intervals <- purrr::map_df(depth.rank.envelope.results, getPinterval)
depth.p.intervals$y.position <- as.numeric(row.names(depth.p.intervals))  # Needed for plotting

# Plot p-value histogram
depth.p.val.hist <- ggplot(data=as.data.frame(depth.p.values), aes(depth.p.values)) + 
  geom_histogram(fill="grey", aes(y = ..ndensity..)) + 
  scale_y_continuous(expand=c(0,0), limits = c(0,1.0)) +
  labs(x = "Rank Envelope p-values", y = "Proportion of trees") +
  theme_procB()
depth.p.val.hist

# Plot p-value ranges
depth.p.int.plot <- ggplot(data = depth.p.intervals, aes(x = Lower, y = y.position, 
                                             xend = Upper, yend=y.position)) 
depth.p.int.plot <- depth.p.int.plot + geom_segment() 
depth.p.int.plot


# Plot MDI statistic histogram #################################################
depth.mdi.val.hist <- ggplot(data=mdi.depth.output, aes(MDI)) + 
  geom_histogram(fill="grey", binwidth = 0.005, aes(y=..ndensity..)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.0)) +
  labs(x = "MDI values", y = "Proportion of Trees") +
  theme_procB()
depth.mdi.val.hist 

# Plot MDI p-values histogram ##################################################
depth.mdi.p.hist <- ggplot(data=mdi.depth.output, aes(MDI.pval)) + 
  geom_histogram(fill="grey", binwidth = 0.001, aes(y=..ndensity..)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.0)) +
  labs(x = "MDI p-values", y = "Proportion of Trees") +
  theme_procB()
depth.mdi.p.hist 
