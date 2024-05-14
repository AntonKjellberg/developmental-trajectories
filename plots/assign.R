library(tidyverse) # Data processing
library(ggpubr) # Plotting and statistics
library(patchwork) # Combine plots
library(speedyseq) # phyloseq with faster functions

# phyloseq object contains raw data from the dada2 pipeline
# Keep in mind that phyloseq naming assumes OTUs. We use ASVs
# otu_table contains ASV counts
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

# Clean and sort time points
sample_data(phy)$Visit <- sample_data(phy)$Visit %>% 
  recode('Week1'='1w', '1mth'='1m', '1yr'='1y', 
         '4yr'  ='4y', '5yr' ='5y', '6yr'='6y') %>% 
  factor(levels = c('1w', '1m', '1y', '4y','5y', '6y'))
sample_data(phy)$Time <- sample_data(phy)$Time %>%
  factor(levels = c('1w', '1m', '1y', '4y', '6y'))

# Filter out rare taxa
tax_table(phy) <- phy %>%  
  transform_sample_counts(function(x) x / sum(x)) %>% 
  filter_taxa(function(x) mean(x) > 0.00001, T)

# Plot presence of each ASV normalized to one
p1 <- phy %>% 
  transform_sample_counts(function(x) x / sum(x)) %>% 
  psmelt() %>% 
  group_by(Species, Time) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  group_by(Species) %>% 
  mutate(Abundance = Abundance / sum(Abundance)) %>% 
  ggplot(aes(y=Abundance, x=Time, group=Species)) +
  geom_line(alpha = 0.3, linewidth = 0.3) +
  theme_bw() +
  theme(legend.position = 'none', axis.title.x=element_blank()) +
  ylab('Standardized ASV Abundance')

# CLR-transform
phy_clr <- phy %>% 
  microbiome::transform('clr')

# Label duplicate ASV names with numbers
tax_table(phy_clr) <- tax_table(phy_clr) %>% 
  data.frame() %>%
  group_by(Species) %>%
  mutate(suffix = row_number()) %>%
  ungroup() %>%
  mutate(Species = ifelse(duplicated(Species), 
                          paste(Species, suffix, sep = ''), Species)) %>%
  select(-suffix) %>% 
  data.frame(row.names = rownames(tax_table(phy_clr))) %>% 
  as.matrix

# Peaks in 1w and 1m:      Assigned 'Early'
# Peaks in 1y:             Assigned 'Intermediate'
# Peaks in 4y and 6y:      Assigned 'Late'
# If 'Late' and abundance at 1y < 2x mean of 4y and 6y: Assigned 'Latest'
colonizers <- phy_clr %>% 
  psmelt() %>%
  group_by(Species, Time) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  group_by(Species) %>%
  mutate(Max_Abundance = max(Abundance)) %>%
  mutate(Colonizer = case_when(
    Time %in% c('1w', '1m')   & Abundance == Max_Abundance ~ 'Early',
    Time  ==    '1y'          & Abundance == Max_Abundance ~ 'Intermediate',
    Time %in% c('4y','6y')    & Abundance == Max_Abundance ~ {
      if (mean(Abundance[Time %in% c('4y','6y')]) >
          2 * Abundance[Time == '1y']) 'Latest' else 'Late' } )) %>% 
  filter(!is.na(Colonizer)) %>% 
  select(Species, Colonizer) %>% 
  distinct(Species, .keep_all = T)

# Add colonizers to the phyloseq object
tax_table(phy_clr) <- tax_table(phy_clr) %>% 
  data.frame() %>% 
  dplyr::left_join(colonizers, by = 'Species') %>%
  data.frame(row.names = rownames(tax_table(phy_clr))) %>% 
  relocate(Colonizer, .before=Domain) %>% 
  as.matrix()

# Make a copy with raw counts and save it. Changes include:
# Trimming out very low abundance ASVs
# Sorting of time points
# Assigning colonizers
tax_table(phy) <- tax_table(phy_clr)
sample_data(phy) <- sample_data(phy_clr)
save(phy, file = 'data/phy_assigned.RData')




# -------------------------------------------------- #
# Plots to characterize the colonizer groups further #
# -------------------------------------------------- #

# Custom legend for the next plot
legend <- tibble(x=1:16, y=1:16, 
                 Colonizer=rep(c('Early','Intermediate', 'Late','Latest'),4)) %>% 
  ggplot(aes(y=y, x=x, group=Colonizer, col=Colonizer)) +
  geom_line(alpha = 0.5, linewidth = 2) +
  scale_color_manual(values=c('darkgreen','orange', 'darkred', 'darkblue')) +
  geom_line(aes(linetype = Colonizer), linewidth = 1, col='black') +
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted', 'dotdash')) +
  theme_bw() +
  theme(legend.key.width = unit(1,'cm'),
        legend.key = element_rect(fill=NA),
        legend.background = element_rect(fill=NA))

# Standardizes values so that they sum to zero
# Retains their standard deviation
standardize0 <- function(x) {(x-mean(x, na.rm=T)) / sum(abs(x), na.rm=T)}

# Plot the standardized colonizer means
plot_means <- tax_glom(phy_clr, taxrank = 'Colonizer') %>%
  psmelt() %>%
  group_by(Colonizer, Time) %>%
  summarize(Mean = mean(Abundance, na.rm = TRUE)) %>%
  group_by(Colonizer) %>% 
  mutate(Mean = standardize0(Mean)) %>% 
  ggplot(aes(y=Mean, x=Time, group=Colonizer, col=Colonizer)) +
  geom_line(linewidth=1) +
  theme_bw() +
  scale_color_manual(values=c('darkgreen','orange', 'darkred', 'darkblue')) +
  ylab('Standardized Colonizer Mean CLR Abundance') +
  theme(legend.position = 'none', axis.title.x=element_blank()) +
  scale_x_discrete(expand = c(0.06,0.06))

# Plot the standardized ASV means
plot_ASVs <- phy_clr %>% 
  psmelt() %>% 
  group_by(Species, Time) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  left_join(tax_table(phy_clr) %>% as_tibble, by = 'Species') %>% 
  group_by(Species) %>% 
  mutate(Abundance = standardize0(Abundance)) %>% 
  ggplot(aes(y=Abundance, x=Time, group=Species, col=Colonizer)) +
  geom_line(linewidth = 0.8, alpha=0.18) +
  theme_bw() +
  scale_color_manual(values=c('darkgreen','orange', 'darkred', 'darkblue')) +
  ylab('Standardized ASV Mean CLR Abundance') +
  stat_summary(aes(group = Colonizer, linetype = Colonizer), fun = mean,
               geom = 'line', linewidth = 1, col = 'black') +
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted', 'dotdash')) +
  theme(legend.position = 'none', axis.title.x=element_blank()) +
  scale_x_discrete(expand = c(0.06,0.06))

# 4x15
plot_means + plot_ASVs + get_legend(legend)

# Clean up species level assignments
# Remove those that are just numbers (hard to interpret)
# Keep only the first two alphabetically (names get too long)
tax_table(phy_clr) <- tax_table(phy_clr) %>% 
  data.frame() %>%
  mutate(Species = gsub('sp[0-9]+/', '', Species),
         Species = gsub('^([^/]+/[^/]+).*', '\\1', Species)) %>%
  data.frame(row.names = rownames(tax_table(phy_clr))) %>% 
  as.matrix

# Plot early ASVs
plot_early <- subset_taxa(phy_clr, Colonizer == 'Early') %>%
  taxa_sums(.) %>% sort %>% tail(12) %>% names %>% 
  prune_taxa(phy_clr) %>% 
  psmelt() %>%
  group_by(Species, Time) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ggplot(aes(y=Abundance, x=Time, group=Species, col=Species, linetype=Species)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  ylab('Mean Relative Log10 Abundance') +
  scale_color_manual(name = 'Early ASVs', 
                     values = rep(c('darkgreen', '#89ab41', '#72C79B'), 4)) +
  scale_linetype_manual(name = 'Early ASVs', 
                        values = rep(c('solid','dashed','dotted','dotdash'), each=3)) +
  theme(legend.title = element_text(face = 'bold'),
        axis.title.x=element_blank(),
        legend.key.width = unit(1,'cm'))

# Plot intermediate ASVs
plot_intermediate <- subset_taxa(phy_clr, Colonizer == 'Intermediate') %>%
  taxa_sums(.) %>% sort %>% tail(12) %>% names %>% 
  prune_taxa(phy_clr) %>%  
  psmelt() %>%
  group_by(Species, Time) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ggplot(aes(y=Abundance, x=Time, group=Species, col=Species, linetype=Species)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  ylab('Mean Relative Log10 Abundance') +
  scale_color_manual(name = 'Intermediate ASVs', 
                     values = rep(c('orange', '#B69006', '#FFDB58'), 4)) +
  scale_linetype_manual(name = 'Intermediate ASVs', 
                        values = rep(c('solid','dashed','dotted','dotdash'), each=3)) +
  theme(legend.title = element_text(face = 'bold'),
        axis.title.x=element_blank(),
        legend.key.width = unit(1,'cm'))

# Plot Late ASVs
plot_late <- subset_taxa(phy_clr, Colonizer == 'Late') %>%
  taxa_sums(.) %>% sort %>% tail(12) %>% names %>% 
  prune_taxa(phy_clr) %>% 
  psmelt() %>%
  group_by(Species, Time) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ggplot(aes(y=Abundance, x=Time, group=Species, col=Species, linetype=Species)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  ylab('Mean Relative Log10 Abundance') +
  scale_color_manual(name = 'Late ASVs', 
                     values = rep(c('darkred', 'red', '#d16680'), 4)) +
  scale_linetype_manual(name = 'Late ASVs', 
                        values = rep(c('solid','dashed','dotted','dotdash'), each=3)) +
  theme(legend.title = element_text(face = 'bold'),
        axis.title.x=element_blank(),
        legend.key.width = unit(1,'cm'))

# Plot Latest ASVs
plot_latest <- subset_taxa(phy_clr, Colonizer == 'Latest') %>%
  taxa_sums(.) %>% sort %>% tail(12) %>% names %>% 
  prune_taxa(phy_clr) %>% 
  psmelt() %>%
  group_by(Species, Time) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ggplot(aes(y=Abundance, x=Time, group=Species, col=Species, linetype=Species)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  ylab('Mean Relative Log10 Abundance') +
  scale_color_manual(name = 'Latest ASVs', 
                     values = rep(c('darkblue', 'cyan4', 'lightblue'), 4)) +
  scale_linetype_manual(name = 'Latest ASVs', 
                        values = rep(c('solid','dashed','dotted','dotdash'), each=3)) +
  theme(legend.title = element_text(face = 'bold'),
        axis.title.x=element_blank(),
        legend.key.width = unit(1,'cm'))

# 15x7
plot_early + plot_intermediate + plot_late + plot_latest + 
  plot_layout(ncol = 1) & theme(legend.justification = 'left')

