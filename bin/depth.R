library('tidyverse')
library('phyloseq')
library('vegan')
library('ggpubr')
library('ggridges')
library('patchwork')
library('phyloseq.extended') # Efficient rarefaction curves

setwd('~/R/copsac/github')

# Phyloseq object contains raw data from the dada2 pipeline
# Keep in mind that phyloseq naming assumes OTUs. We use ASVs
# otu_matrix contains ASV counts
load('data/dada2_gtdb_2021_11_18.RData')

# Clean up and Sort time points
sample_data(phy)$Visit <- sample_data(phy)$Visit %>% 
  recode('Week1'='1w', '1mth'='1m', '1yr'='1y', 
         '4yr'  ='4y', '5yr' ='5y', '6yr'='6y') %>% 
  factor(levels = c('1w', '1m', '1y', '4y','5y', '6y'))

# Only include fecal samples from the child
# Remove taxa that are not present
phy <- phy %>% 
  subset_samples(Visit %in% c('1w', '1m', '1y', '4y', '5y', '6y')) %>%
  subset_samples(Sampletype == 'Faeces') %>% 
  prune_taxa(taxa_sums(.) > 0, .)

# Minimum and mean depths
sample_data(phy)$Reads %>% mean %>% round
sample_data(phy)$Reads %>% min

# Colors for time points
cols <- c('#f65151', '#FFA351', '#2CA674', '#1875B9', '#7f60ab', '#D079B1')

# Legend for the following plots
legend <- data.frame(Timepoint = unique(sample_data(phy)$Visit), Values = 1:6,
                     Lines = rep(c('Trim: 2020', 'Mean: 54748       '), each=3)) %>% 
  ggplot(aes(x=Values, y=Values)) +
  geom_point(aes(fill=Timepoint), shape = 21, alpha = 0.7, size=2) +
  geom_line(aes(color=Lines), linetype='dashed') +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = c('red', 'black')) +
  theme_minimal() +
  theme(legend.position='top',
        legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 1))

# Depth distribution of all samples
plot_depth_all <- sample_data(phy) %>% 
  ggplot(aes(x=Reads)) +
  geom_histogram(bins = 50, fill = 'lightgrey', alpha = 0.7) +
  geom_vline(xintercept = mean(sample_data(phy)$Reads), 
             linetype = 'dashed', col = 'red') +
  ylab('Frequency') +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 100000, 200000), 
                     labels = c(0, '100k', '200k'),
                     limits = c(0, 200000))

# Depth distributions grouped by time point
plot_depth_time <- sample_data(phy) %>% 
  ggplot(aes(x=Reads, y=Visit, fill=Visit)) +
  geom_density_ridges(aes(point_fill=Visit, point_shape=21),
                      point_alpha = 0.1, alpha = 0.5, jittered_points = T) +
  scale_discrete_manual(aesthetics = c('point_fill', 'fill'), values = cols) + 
  scale_x_continuous(breaks = c(0, 100000, 200000), 
                     labels = c(0, '100k', '200k'),
                     limits = c(0, 200000)) +
  theme_bw() +
  theme(legend.position = 'none', 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# Generate rarefaction curves for all samples
# This takes ~20m on a good laptop
rarefaction_curves <- ggrare(phy, step = 200, se = FALSE, color = 'Visit')

# Rarefaction curves of all samples
plot_rare_all <- rarefaction_curves$data %>%
  ggplot(aes(x=Size, y = .S, group = Sample, col = Visit)) +
  geom_line(alpha = 0.3, linewidth = 0.05) +
  geom_vline(xintercept = 2020, linetype = 'dashed') +
  scale_x_continuous(breaks = c(0, 10000, 20000), 
                     labels = c(0, '10k', '20k'),
                     limits = c(0, 20000)) +
  scale_color_manual(values = cols) +
  theme_bw() +
  xlab('Sequences per sample') +
  ylab('ASV Richness') +
  theme(legend.position = 'none')

# rarefaction curves grouped by time point
plot_rare_separate <- rarefaction_curves$data %>%
  mutate(Size = log(Size)) %>% 
  ggplot(aes(x=Size, y = .S, group = Sample, col = Visit)) +
  facet_wrap(~Visit) +
  geom_line(alpha = 0.5, linewidth = 0.65) +
  scale_color_manual(values = cols) +
  geom_vline(xintercept = log(2020), linetype = 'dashed') +
  theme_bw() +
  xlab('Log10 sequences per sample') +
  ylab('ASV Richness') +
  theme(legend.position = 'none', strip.text.x = element_blank())

# 1000x1000
(as_ggplot(get_legend(legend)) / 
    ((((plot_depth_all / plot_rare_all) | plot_depth_time) +
        plot_layout(widths = c(2,1)))) / plot_rare_separate + 
    plot_layout(nrow = 3)) + plot_layout(heights = c(1, 10, 10))


# Extract relevant metadata for samples
metadata <- phy %>% 
  sample_data() %>% 
  as_tibble() %>% 
  select(SampleID, Reads, Visit)

# Generate a otu table in long format
otu_long <- phy %>%
  otu_table() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'ASV') %>% 
  pivot_longer(-ASV, names_to='Sample', values_to='Raw')

# Add column with rarefied measures (Vegan is a lot faster than phyloseq)
# Combine with relevant metadata
# Calculate alpha diversity measures 
otu_long <- phy %>% 
  otu_table() %>% 
  as.data.frame() %>% 
  t() %>% rrarefy(min(metadata$Reads)) %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'ASV') %>% 
  pivot_longer(-ASV, names_to='Sample', values_to='Subsampled') %>% 
  mutate(Raw = otu_long$Raw) %>% 
  pivot_longer(cols = c(Raw, Subsampled), names_to = 'Subsampled', 
               values_to = 'Reads') %>% 
  group_by(Sample, Subsampled) %>% 
  summarize(Shannon = diversity(Reads, index='shannon'),
            Richness = sum(Reads > 0)) %>% 
  dplyr::left_join(metadata, by=c('Sample' = 'SampleID')) %>% 
  pivot_longer(cols = c(Shannon, Richness), names_to = 'Diversity', 
               values_to = 'Value')

# Richness diversity by reads
plot_richness <- otu_long %>% 
  filter(Diversity == 'Richness') %>% 
  ggscatter(x='Reads', y='Value', col='Visit', size = 0.85,
            add='reg.line', alpha=0.1, conf.int=T, palette=cols) +
  facet_wrap(~Subsampled) +
  xlim(0, 200000) +
  ylim(0, 400) +
  theme_bw() +
  stat_cor(aes(color = Visit), method = 'spearman', size = 3.5,
           label.x.npc = 'left', label.y.npc = 'top') +
  ylab('ASV Richness') +
  theme(legend.position = 'none',
        strip.text = element_text(size = 11),
        strip.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Shannon diversity by reads
plot_shannon <- otu_long %>% 
  filter(Diversity == 'Shannon') %>% 
  ggscatter(x='Reads', y='Value', col='Visit', fullrange = F, size = 0.85, 
            add='reg.line', alpha=0.1, conf.int=T, palette=cols) +
  facet_wrap(~Subsampled) +
  xlim(0, 200000) +
  ylim(0, 8) +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 50000, 100000, 150000, 200000), 
                     labels = c(0, '50k', '100k', '150k', '200k'),
                     limits = c(0, 200000)) +
  stat_cor(aes(color = Visit), method = 'spearman', size = 3.5,
           label.x.npc = 'left', label.y.npc = 'top') +
  ylab('Shannon Index') +
  theme(legend.position = 'none',
        strip.text.x = element_blank())

# 800x700
plot_richness / plot_shannon

