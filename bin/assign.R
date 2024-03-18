library('tidyverse')
library('phyloseq') # I use speedyseq for slow phyloseq functions
library('ggpubr')
library('patchwork')

setwd('~/R/copsac/github')

# Phyloseq object contains raw data from the dada2 pipeline
# Keep in mind that phyloseq naming assumes OTUs. We use ASVs
# otu_matrix contains ASV counts
load('data/dada2_gtdb_2021_11_18.RData')

# Exclude irrelevant samples
# Visit: keeps 5y time point
# Time: allocates 5y samples to 4y or 6y
# 4y if no sample at 4y already. Otherwise 6y
sample_data(phy) <- phy %>% 
  subset_samples(Time %in% c('1w', '1m', '1y', '4y', '5y', '6y')) %>%
  subset_samples(Sampletype == 'Faeces') %>% 
  sample_data %>% 
  data.frame %>% 
  mutate(Time = case_when(
    Time == '5y' & ABCNO %in% filter(., Time == '4y')$ABCNO ~ '6y',
    Time == '5y' ~ '4y',
    TRUE ~ Time))

# Transform to relative abundance
# Filter out very low abundance ASVs
phy <- phy %>%  
  microbiome::transform('compositional') %>% 
  filter_taxa(function(x) mean(x) > 0.0000001, T)

# Clean and Sort time points
sample_data(phy)$Visit <- sample_data(phy)$Visit %>% 
  recode('Week1'='1w', '1mth'='1m', '1yr'='1y', 
         '4yr'  ='4y', '5yr' ='5y', '6yr'='6y') %>% 
  factor(levels = c('1w', '1m', '1y', '4y','5y', '6y'))
sample_data(phy)$Time <- sample_data(phy)$Time %>%
  factor(levels = c('1w', '1m', '1y', '4y', '6y'))

# Label duplicate ASV names with numbers to differentiate them
tax_table(phy) <- tax_table(phy) %>% 
  data.frame() %>%
  group_by(Species) %>%
  mutate(suffix = row_number()) %>%
  ungroup() %>%
  mutate(Species = ifelse(duplicated(Species), 
                          paste(Species, suffix, sep = ''), Species)) %>%
  select(-suffix) %>% 
  data.frame(row.names = rownames(tax_table(phy))) %>% 
  as.matrix

# Assign colonizer groups to the ASVs
# Peaks in 1w and 1m:      Assigned 'Early'
# Peaks in 1y:             Assigned 'Intermediate'
# Peaks in 4y, 5y, and 6y: Assigned 'Late'
# If 'Late' and abundance at 1y < 3x mean of 4y, 5y, and 6y: Assigned 'Latest'
colonizers <- phy %>% 
  speedyseq::psmelt() %>%
  group_by(Species, Visit) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  group_by(Species) %>%
  mutate(Max_Abundance = max(Abundance)) %>%
  mutate(Colonizer = case_when(
    Visit %in% c('1w', '1m')   & Abundance == Max_Abundance ~ 'Early',
    Visit  ==    '1y'          & Abundance == Max_Abundance ~ 'Intermediate',
    Visit %in% c('4y','5y','6y')    & Abundance == Max_Abundance ~ {
      if (mean(Abundance[Visit %in% c('4y','5y','6y')]) >
          3 * Abundance[Visit == '1y']) 'Latest' else 'Late' } )) %>% 
  filter(!is.na(Colonizer)) %>% 
  select(Species, Colonizer) %>% 
  distinct(Species, .keep_all = T)

# Add colonizers to the phyloseq object
tax_table(phy) <- tax_table(phy) %>% 
  data.frame() %>% 
  dplyr::left_join(colonizers, by = 'Species') %>%
  data.frame(row.names = rownames(tax_table(phy))) %>% 
  relocate(Colonizer, .before=Domain) %>% 
  as.matrix()

# Save the phyloseq. Changes include:
# Mean Relative abundance transformation
# Trimming out very low abundance ASVs
# Sorting of time points
# Assigning colonizers
save(phy, file = 'data/phy_assigned.RData')

# Plot the Colonizer means
tax_glom(phy, taxrank = "Colonizer") %>%
  speedyseq::psmelt() %>%
  group_by(Colonizer, Visit) %>%
  summarize(Mean = mean(Abundance, na.rm = TRUE)) %>%
  ggplot(aes(y=Mean, x=Visit, group=Colonizer, col=Colonizer)) +
  geom_line(linewidth=1) +
  theme_bw() +
  scale_color_manual(values=c('darkgreen','orange', 'darkred', 'darkblue')) +
  ylab('Mean Relative Abundance')

# Custom legend for the next two plots
legend <- phy %>% 
  filter_taxa(function(x) mean(x) > 0.0005, T) %>% 
  speedyseq::psmelt() %>% 
  group_by(Species, Visit) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  left_join(colonizers, by = 'Species') %>% 
  ggplot(aes(y=Abundance, x=Visit, group=Species, col=Colonizer)) +
  geom_line(alpha = 0.5, linewidth = 2) +
  theme_bw() +
  scale_color_manual(values=c('darkgreen','orange', 'darkred', 'darkblue')) +
  stat_summary(aes(group = Colonizer, linetype = Colonizer), fun = mean,
               geom = 'line', linewidth = 1, col = 'black') +
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted', 'dotdash')) +
  theme(legend.key.width = unit(1,'cm'))

# Transformed ASV abundance means
# Used to fade out low abundance ASVs in the next plot
alpha <- phy %>% 
  filter_taxa(function(x) mean(x) > 0.0005, T) %>% 
  speedyseq::psmelt() %>% 
  group_by(Species) %>% 
  summarize(Abundance = sum(Abundance)) %>% 
  mutate(Abundance = Abundance / max(Abundance)) %>% 
  .[[2]] %>% rep(each = 6) 

# Plot presence of each ASV normalized to one
plot_normalized <- phy %>% 
  filter_taxa(function(x) mean(x) > 0.0005, T) %>% 
  speedyseq::psmelt() %>% 
  group_by(Species, Visit) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  left_join(colonizers, by = 'Species') %>% 
  group_by(Species) %>% 
  mutate(Abundance = Abundance / sum(Abundance)) %>% 
  ggplot(aes(y=Abundance, x=Visit, group=Species, col=Colonizer)) +
  geom_line(alpha = sqrt(alpha)*1.5, linewidth = 0.8) +
  theme_bw() +
  scale_color_manual(values=c('darkgreen','orange', 'darkred', 'darkblue')) +
  ylab('Normalized ASV Abundance') +
  stat_summary(aes(group = Colonizer, linetype = Colonizer), fun = mean,
               geom = 'line', linewidth = 1.5, col = 'black') +
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted', 'dotdash')) +
  theme(legend.position = 'none', axis.title.x=element_blank())

plot_normalized + plot_spacer() + get_legend(legend) + 
  plot_layout(widths = c(5, -1.7, 5))

# Z-score normalization function
z_score <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Plot z-score normalized abundance of each ASV
plot_z_score <- phy %>% 
  speedyseq::psmelt() %>% 
  mutate(Visit = case_when(
    Visit %in% c('1w', '1m') ~ 'Before',
    Visit  ==    '1y'        ~ '1 year',
    TRUE                     ~ 'After')) %>% 
  mutate(Visit = Visit %>% factor(levels = c('Before', '1 year', 'After'))) %>% 
  group_by(Species, Visit) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  left_join(colonizers, by = 'Species') %>% 
  group_by(Species) %>% 
  mutate(Abundance = z_score(Abundance)) %>% 
  ggplot(aes(y=Abundance, x=Visit, group=Species, col=Colonizer)) +
  geom_line(linewidth = 0.7, alpha=0.2) +
  theme_bw() +
  scale_color_manual(values=c('darkgreen','orange', 'darkred', 'darkblue')) +
  ylab('ASV Abundance Z-score') +
  stat_summary(aes(group = Colonizer, linetype = Colonizer), fun = mean,
               geom = 'line', linewidth = 1, col = 'black') +
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted', 'dotdash')) +
  theme(legend.position = 'none', axis.title.x=element_blank())

plot_z_score + plot_spacer() + get_legend(legend) + 
  plot_layout(widths = c(5, -1.7, 5))

# Clean up species level assignments
# Remove those that are just numbers (hard to interpret)
# Keep only the first two alphabetically (names get too long)
tax_table(phy) <- tax_table(phy) %>% 
  data.frame() %>%
  mutate(Species = gsub('sp[0-9]+/', '', Species),
         Species = gsub('^([^/]+/[^/]+).*', '\\1', Species)) %>%
  data.frame(row.names = rownames(tax_table(phy))) %>% 
  as.matrix

# Plot early ASVs
plot_early <- subset_taxa(phy, Colonizer == 'Early') %>%
  filter_taxa(function(x) mean(x) > 0.007, T) %>% 
  psmelt() %>%
  group_by(Species, Visit) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ggplot(aes(y=Abundance, x=Visit, group=Species, col=Species, linetype=Species)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  ylab('Mean Relative Abundance') +
  scale_color_manual(name = 'Early ASVs', 
                     values = rep(c('darkgreen', '#89ab41', '#72C79B'), 4)) +
  scale_linetype_manual(name = 'Early ASVs', 
                        values = rep(c('solid','dashed','dotted','dotdash'), each=3)) +
  theme(legend.title = element_text(face = 'bold'),
        axis.title.x=element_blank())

# Plot intermediate ASVs
plot_intermediate <- subset_taxa(phy, Colonizer == 'Intermediate') %>%
  filter_taxa(function(x) mean(x) > 0.004, T) %>% 
  psmelt() %>%
  group_by(Species, Visit) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ggplot(aes(y=Abundance, x=Visit, group=Species, col=Species, linetype=Species)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  ylab('Mean Relative Abundance') +
  scale_color_manual(name = 'Intermediate ASVs', 
                     values = rep(c('orange', '#B69006', '#FFDB58'), 4)) +
  scale_linetype_manual(name = 'Intermediate ASVs', 
                        values = rep(c('solid','dashed','dotted','dotdash'), each=3)) +
  theme(legend.title = element_text(face = 'bold'),
        axis.title.x=element_blank())

# Plot Late ASVs
plot_late <- subset_taxa(phy, Colonizer == 'Late') %>%
  filter_taxa(function(x) mean(x) > 0.004, T) %>% 
  psmelt() %>%
  group_by(Species, Visit) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ggplot(aes(y=Abundance, x=Visit, group=Species, col=Species, linetype=Species)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  ylab('Mean Relative Abundance') +
  scale_color_manual(name = 'Late ASVs', 
                     values = rep(c('darkred', 'red', '#d16680'), 4)) +
  scale_linetype_manual(name = 'Late ASVs', 
                        values = rep(c('solid','dashed','dotted','dotdash'), each=3)) +
  theme(legend.title = element_text(face = 'bold'),
        axis.title.x=element_blank())

# Plot Latest ASVs
plot_latest <- subset_taxa(phy, Colonizer == 'Latest') %>%
  filter_taxa(function(x) mean(x) > 0.004, T) %>% 
  psmelt() %>%
  group_by(Species, Visit) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ggplot(aes(y=Abundance, x=Visit, group=Species, col=Species, linetype=Species)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  ylab('Mean Relative Abundance') +
  scale_color_manual(name = 'Latest ASVs', 
                     values = rep(c('darkblue', 'cyan4', 'lightblue'), 4)) +
  scale_linetype_manual(name = 'Latest ASVs', 
                        values = rep(c('solid','dashed','dotted','dotdash'), each=3)) +
  theme(legend.title = element_text(face = 'bold'),
        axis.title.x=element_blank())

plot_early + plot_intermediate + plot_late + plot_latest + 
plot_layout(ncol = 1) & theme(legend.justification = 'left')

