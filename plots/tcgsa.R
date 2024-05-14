library(tidyverse) # Data processing
library(speedyseq) # Phyloseq with faster functions
library(TcGSA) # Time-course gene set analysis
library(patchwork) # Combine plots

load('data/dada2_gtdb_2021_11_18.RData')

# Split 5y into 4y and 6y. 4y if no sample at 4y already. Otherwise 5y
sample_data(phy) <- phy %>% 
  sample_data %>% 
  data.frame %>% 
  mutate(Time = case_when(
    Time == '5y' & ABCNO %in% filter(., Time == '4y')$ABCNO ~ '6y',
    Time == '5y' ~ '4y',
    TRUE ~ Time))

# Clean up time points
sample_data(phy)$Visit <- sample_data(phy)$Visit %>% 
  recode('Week1'='1w', '1mth'='1m', '1yr'='1y', 
         '4yr'  ='4y', '5yr' ='5y', '6yr'='6y') %>% 
  factor(levels = c('1w', '1m', '1y', '4y', '5y', '6y'))
sample_data(phy)$Time <- sample_data(phy)$Time %>% 
  factor(levels = c('1w', '1m', '1y', '4y', '6y'))

# Only include fecal samples from the child
phy <- phy %>% 
  subset_samples(Visit %in% c('1w', '1m', '1y', '4y', '5y', '6y')) %>%
  subset_samples(Sampletype == 'Faeces')




# -------------------------------------------------------------- #
# Assigning colonizer groups using Time-Course Gene Set Analysis #
# -------------------------------------------------------------- #


# Agglomerate to genus level and transform to relative abundance
phy <- tax_glom(phy, 'Genus') %>% 
  transform_sample_counts(function(x) x / sum(x))


# Subset of genera present in more than 20% of samples
# With a mean relative abundance above 0.001% of samples at a time point
genera <- c()
for (visit in levels(sample_data(phy)$Visit)) {
  subset <- phy %>% 
    subset_samples(Visit == visit) %>% 
    filter_taxa(function(x) sum(x > 0) / length(x) >= 0.2, TRUE)
  # Calculate mean relative abundance
  subset <- data.frame(Genus = as.character(tax_table(subset)[,'Genus']),
                       mra = taxa_sums(subset) / nsamples(subset))
  genera <- c(genera, filter(subset, mra > 1e-3) %>% pull(Genus)) %>% unique
}
phy <- prune_taxa(tax_table(phy)[,'Genus'] %in% genera, phy)

# Log-transform ratios. Pseudocount is mean relative abundance cutoff
# I've tested combinations of Z-score, CLR, Log10, MRA, and ranks
phy <- transform_sample_counts(phy, function(x) log(1e-3+x))

# Extract data for the TcGSA clustering
otu_table <- phy %>% 
  otu_table() %>% 
  data.frame(row.names = as_tibble(tax_table(phy))$Genus)
genera <- list(genesets = list(as_tibble(tax_table(phy))$Genus),
               geneset.names = 'AllGenera',
               geneset.descriptions = 'All genera')
IDs <- sample_data(phy) %>% row.names()

# I've tried scaling the time points differently
# 1w and 1m time points are swapped in the original script
# I decided to correct it here. No influence on the clustering
time <- sample_data(phy)$Visit %>% recode('1w'=1,'1m'=2,'1y'=3,
                                          '4y'=4,'5y'=5, '6y'=6)

# Use TcGSA to assign colonizers
# aggreg.fun = 'mean' is not used by Roswall et al
trajectories <- plot1GS(expr=otu_table, gmt=genera,
                        Subject_ID = IDs, TimePoint = time,
                        geneset.name='AllGenera',
                        aggreg.fun = 'mean')

# Add trajectories to the phyloseq taxonomy table
tax_table(phy) <- tax_table(phy) %>% 
  data.frame() %>% 
  left_join(trajectories$classif, by = c('Genus' = 'ProbeID')) %>%
  data.frame(row.names = rownames(tax_table(phy))) %>% 
  rename(Colonizer = 'Cluster') %>% 
  relocate(Colonizer, .before=Domain) %>% 
  as.matrix() %>%
  tax_table

# legend for the next plot
legend <- tibble(x=1:16, y=1:16, Trajectory=factor(rep(c(1:4),4))) %>% 
  ggplot(aes(y=y, x=x, group=Trajectory, col=Trajectory)) +
  geom_line(alpha = 0.5, linewidth = 2) +
  scale_color_manual(values=c('darkgreen','orange', 'darkred', 'darkblue')) +
  geom_line(aes(linetype = Trajectory), linewidth = 1, col='black') +
  scale_linetype_manual(values = c('solid', 'dotdash', 'dashed', 'dotted')) +
  guides(color = guide_legend(nrow = 2)) +
  theme(legend.key.width = unit(1,'cm'),
        legend.key = element_rect(fill=NA),
        legend.background = element_rect(fill=NA))

# Plot the TcGSA estimates (mixed effects model)
estimates <- trajectories$p$data %>%
  mutate(Time = case_when(
    TimePoint == 1 ~ '1w',
    TimePoint == 2 ~ '1m',
    TimePoint == 3 ~ '1y',
    TimePoint == 4 ~ '4y',
    TimePoint == 5 ~ '5y',
    TimePoint == 6 ~ '6y')) %>%
  mutate(Time = factor(Time, levels = c('1w', '1m', '1y', '4y', '5y', '6y'))) %>%
  ggplot(aes(x=Time, y=value, col=Cluster, group=Probe_ID)) +
  geom_line(alpha=0.4, linewidth = 0.7) +
  stat_summary(aes(linetype=Cluster, group=Cluster),
               fun='median', geom='line', col='black', linewidth=0.9) +
  scale_linetype_manual(values = c('solid', 'dotted', 'dotdash', 'dashed')) +
  scale_color_manual(values=c('darkgreen','darkblue','orange','darkred')) +
  theme_bw() +
  ylab('Mean of Standardized Estimate') +
  theme(legend.position = 'none', axis.title.x=element_blank())

estimates + inset_element(as_ggplot(get_legend(legend)),0.7,0,0,1.6)




# ---------------------------------- #
# Plots to characterize trajectories #
# ---------------------------------- #

# I swap the names of the colonizers to make it more easy to interpret
# (1 = early, 2 = intermediate, 3 = late, 4 = latest)
plot1 <- subset_taxa(phy, Colonizer == 1) %>%
  taxa_sums(.) %>% sort %>% tail(12) %>% names %>% 
  prune_taxa(phy) %>%  
  psmelt() %>%
  group_by(Genus, Visit) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ggplot(aes(y=Abundance, x=Visit, group=Genus, col=Genus, linetype=Genus)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  ylab('Mean Log10 Relative Abundance') +
  scale_color_manual(name = 'Trajectory 1 Genera', 
                     values = rep(c('darkgreen', '#89ab41', '#72C79B'), 4)) +
  scale_linetype_manual(name = 'Trajectory 1 Genera', 
                        values = rep(c('solid','dashed','dotted','dotdash'), each=3)) +
  theme(legend.title = element_text(face = 'bold'),
        axis.title.x=element_blank(),
        legend.key.width = unit(1,'cm'))

plot2 <- subset_taxa(phy, Colonizer == 3) %>%
  taxa_sums(.) %>% sort %>% tail(12) %>% names %>% 
  prune_taxa(phy) %>% 
  psmelt() %>%
  group_by(Genus, Visit) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ggplot(aes(y=Abundance, x=Visit, group=Genus, col=Genus, linetype=Genus)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  ylab('Mean Log10 Relative Abundance') +
  scale_color_manual(name = 'Trajectory 2 Genera', 
                     values = rep(c('orange', '#B69006', '#FFDB58'), 4)) +
  scale_linetype_manual(name = 'Trajectory 2 Genera', 
                        values = rep(c('solid','dashed','dotted','dotdash'), each=3)) +
  theme(legend.title = element_text(face = 'bold'),
        axis.title.x=element_blank(),
        legend.key.width = unit(1,'cm'))

plot3 <- subset_taxa(phy, Colonizer == 4) %>%
  taxa_sums(.) %>% sort %>% tail(12) %>% names %>% 
  prune_taxa(phy) %>% 
  psmelt() %>%
  group_by(Genus, Visit) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ggplot(aes(y=Abundance, x=Visit, group=Genus, col=Genus, linetype=Genus)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  ylab('Mean Log10 Relative Abundance') +
  scale_color_manual(name = 'Trajectory 3 Genera', 
                     values = rep(c('darkred', 'red', '#d16680'), 4)) +
  scale_linetype_manual(name = 'Trajectory 3 Genera', 
                        values = rep(c('solid','dashed','dotted','dotdash'), each=3)) +
  theme(legend.title = element_text(face = 'bold'),
        axis.title.x=element_blank(),
        legend.key.width = unit(1,'cm'))

plot4 <- subset_taxa(phy, Colonizer == 2) %>%
  taxa_sums(.) %>% sort %>% tail(12) %>% names %>% 
  prune_taxa(phy) %>% 
  psmelt() %>%
  group_by(Genus, Visit) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ggplot(aes(y=Abundance, x=Visit, group=Genus, col=Genus, linetype=Genus)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  ylab('Mean Log10 Relative Abundance') +
  scale_color_manual(name = 'Trajectory 4 Genera', 
                     values = rep(c('darkblue', 'cyan4', 'lightblue'), 4)) +
  scale_linetype_manual(name = 'Trajectory 4 Genera', 
                        values = rep(c('solid','dashed','dotted','dotdash'), each=3)) +
  theme(legend.title = element_text(face = 'bold'),
        axis.title.x=element_blank(),
        legend.key.width = unit(1,'cm'))

# 15x7
plot1 + plot2 + plot3 + plot4 + 
  plot_layout(ncol = 1) & theme(legend.justification = 'left')
