# 08_variance_analysis.R

# Author:  Kara Feilich
# Date modified: 4 September 2018

# Looking at the variance in PC scores accrued before and after a narrow Eocene 
# time window.

# Find variance in PC scores for earliest fossils
in_window_pc_variance <- pelagia.pcscores.wHTs %>%
  filter(InTimeWindow == 1) %>%
  dplyr::select(starts_with("PC")) %>%
  summarise_all(var) %>%
  rowSums()

# Find variance in PC scores for fossils that aren't in the window
not_window_pc_variance <- pelagia.pcscores.wHTs %>%
  filter(InTimeWindow == 0 & Status == "fossil") %>%
  dplyr::select(starts_with("PC")) %>%
  summarise_all(var) %>%
  rowSums()

modern_pc_variance <- pelagia.pcscores.wHTs %>%
  filter(Status == "modern") %>%
  dplyr::select(starts_with("PC")) %>%
  summarise_all(var) %>%
  rowSums()

fossil_pc_variance <- pelagia.pcscores.wHTs %>%
  filter(Status == "fossil") %>%
  dplyr::select(starts_with("PC")) %>%
  summarise_all(var) %>%
  rowSums()


print(in_window_pc_variance/modern_pc_variance)
print(not_window_pc_variance/modern_pc_variance)
print(fossil_pc_variance/modern_pc_variance)
