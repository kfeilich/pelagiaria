# 03_plot_DTT_figures.R

# Author: Kara Feilich
# Date modified: 08 June 2018
# Run this after 02_DTT_analysis.R to plot interesting things
load(file = "data/AR.2tailed.rankenv.RData")


# Derive MDI values ##########################################################
# Pull the values of the MDI statistic into a vector
mdi.values <- purrr::map_dbl(dtt.results.all, getMDI)
# Pull the p-values associated with the MDI statistic into a vector
mdi.p.values <- purrr::map_dbl(dtt.results.all, getMDIPvalue)
# Bind the MDI values and the p-values into a dataframe
mdi.output <- data.frame("MDI" = mdi.values, "MDI.pval" = mdi.p.values)


# Make rank envelope p-value range histograms ##################################
# Pulls all p-values from list of rank-envelope results into a vector (type double)
p.values <- purrr::map_dbl(AR.envelope.results, getPvalue)

# Pulls all p-value intervals from list of rank envelope results into a dataframe
# with columns of Lower (bound) and Upper (bound)
p.intervals <- purrr::map_df(AR.envelope.results, getPinterval)
p.intervals$y.position <- as.numeric(row.names(p.intervals))  # Needed for plotting

# Plot p-value histogram
p.val.hist <- ggplot(data=as.data.frame(p.values), aes(p.values)) + 
  geom_histogram(fill="grey", binwidth = 0.0025, aes(y = ..ndensity..)) + 
  scale_y_continuous(expand=c(0,0), limits = c(0,1.0)) +
  labs(x = "Rank Envelope p-values", y = "Proportion of trees") +
  theme_procB()
p.val.hist

# Plot p-value ranges
p.int.plot <- ggplot(data = p.intervals, aes(x = Lower, y = y.position, 
                                             xend = Upper, yend=y.position)) 
p.int.plot <- p.int.plot + geom_segment() 
p.int.plot


# Plot MDI statistic histogram #################################################
mdi.val.hist <- ggplot(data=mdi.output, aes(MDI)) + 
  geom_histogram(fill="grey", binwidth = 0.005, aes(y=..ndensity..)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.0)) +
  labs(x = "MDI values", y = "Proportion of Trees") +
  theme_procB()
mdi.val.hist 

# Plot MDI p-values histogram ##################################################
mdi.p.hist <- ggplot(data=mdi.output, aes(MDI.pval)) + 
  geom_histogram(fill="grey", binwidth = 0.001, aes(y=..ndensity..)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.0)) +
  labs(x = "MDI p-values", y = "Proportion of Trees") +
  theme_procB()
mdi.p.hist 

rm(dtt.results.all, AR.envelope.results, cleaned.trees, agg.timepts, sig.bins,
   mdi.output, mdi.p.hist, mdi.val.hist, p.int.plot, p.intervals, p.val.hist)