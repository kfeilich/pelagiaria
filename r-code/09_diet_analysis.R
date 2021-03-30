# 02b_diet_analysis.R

# Author: Kara Feilich
# Date modified: 19 February 2019

trees <- read.nexus ("data/Posterior_100_pruned.nex")
diet_colors <- c("fish" = "black", "jellyfish" = "dodgerblue", "zooplankton" = "green")

#Import and prepare trait data for analysis
diet <- readr::read_csv("data/Pelagia_diet_data.csv") %>%
  dplyr::select(Family, Species, "Zoo-Jelly-Fish") %>%
  dplyr::rename(Zoo_Jelly_Fish = "Zoo-Jelly-Fish") %>%
  dplyr::mutate(Zoo_Jelly_Fish = as.factor(Zoo_Jelly_Fish),
         Species = str_to_title(str_extract(string = Species, pattern = "^[:alpha:]+_[:alpha:]+"))) %>%
  na.omit
diet
  
diet_trait <- pull(diet,Zoo_Jelly_Fish) 
names(diet_trait) <- dplyr::pull(diet, Species)  

set.seed(7)

# Runs CTT on list of trees
ctt_outputs <- map2(trees,c(1:length(trees)), do_ctt, diet_trait, n_sim=100)
save(ctt_outputs, file="data/ctt_outputs.RData")

# Compares observed values to null envelope
ctt_envelopes_rate <- map(ctt_outputs, check_ctt_env, type="rate")
save(ctt_envelopes_rate, file = "data/ctt_envelopes_rates.RData")

# Returns vectors indicating if data is above or below envelope
ctt_envelopes_number <- map(ctt_outputs, check_ctt_env, type="number")
save(ctt_envelopes_number, file = "data/ctt_envelopes_number.RData")

rm(ctt_outputs) # Remove from env. to save RAM

ctt_sig_times_number <- map(ctt_envelopes_number, lower_or_higher) # Find number lower or higher
ctt_sig_times_rate <- map(ctt_envelopes_rate, lower_or_higher) # Find number lower or higher

#Get the number of bins in each analysis
ctt_sig_times_nbins <- map_dbl(ctt_sig_times_number,length)
max_ctt_nbins <- max(ctt_sig_times_nbins)  # Max number of bins
analyses_to_0pad <- which(ctt_sig_times_nbins != max_ctt_nbins)
names(analyses_to_0pad) <- NULL

rm(ctt_envelopes_number) # Remove from env. to save RAM

# load("data/ctt_envelopes_rates.RData")
ctt_sig_times_rates <- map(ctt_envelopes_rate, lower_or_higher)
rm(ctt_envelopes_rate) # Remove from env. to save RAM


# Extend analyses with bins of 0 up to max bins
for (i in analyses_to_0pad) { 
  ctt_sig_times_number[[i]] <- c(ctt_sig_times_number[[i]], 
                                 replicate(max_ctt_nbins - length(ctt_sig_times_number[[i]]), 0))
  ctt_sig_times_rates[[i]] <- c(ctt_sig_times_rates[[i]], 
                                 replicate(max_ctt_nbins - length(ctt_sig_times_rates[[i]]), 0))
  }
  


number_by_timepoint <- t(bind_cols(ctt_sig_times_number))
rate_by_timepoint <- t(bind_cols(ctt_sig_times_rates))

number_above_timepoint <- data.frame(number_by_timepoint) %>%
  summarize_all( function(x) length(which(x == 1)))
number_below_timepoint <- data.frame(number_by_timepoint) %>%
  summarize_all( function(x) length(which(x == -1)))

number_changes_plotting <- data.frame(Bins = c(1:max_ctt_nbins), 
                                      number_above = as.numeric(number_above_timepoint),
                                      number_below = as.numeric(number_below_timepoint))

rate_above_timepoint <- data.frame(rate_by_timepoint) %>%
  summarize_all( function(x) length(which(x == 1)))
rate_below_timepoint <- data.frame(rate_by_timepoint) %>%
  summarize_all( function(x) length(which(x == -1)))

rate_changes_plotting <- data.frame(Bins = c(1:max_ctt_nbins), 
                                      rate_above = as.numeric(rate_above_timepoint),
                                      rate_below = as.numeric(rate_below_timepoint))

number_changes_plot <- ggplot(data = number_changes_plotting, aes(x=Bins)) + 
  geom_bar(aes(y=number_above), stat =  "identity", fill = "black") +
  geom_bar(aes(y=number_below), stat =  "identity", fill = "black") +
  xlab("Time Bin (million years ago)") +
  ylab("Number Posterior Trees") +
  ylim(c(-100,100))+
  scale_x_reverse(limits=c(max_ctt_nbins,0), expand = c(0,0))+
  ggtitle("Diet")+
  theme_custom()