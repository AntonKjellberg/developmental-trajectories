library(tidyverse) # Data processing
library(speedyseq) # Phyloseq with faster functions
library(vegan) # Community ecology processing
library(ranacapa) # Fast rarefaction curves
library(ggpubr) # Statistics and tables in plots
library(gghalves) # Half-half plots
library(ggridges) # Ridgeline plots
library(patchwork) # Combine plots
library(rstatix)  # Pipe-friendly statistic tests

# Phyloseq object contains raw data from the dada2 pipeline
# Keep in mind that phyloseq naming assumes OTUs. We use ASVs
# otu_matrix contains ASV counts
load('data/dada2_gtdb_2021_11_18.RData')

# Clean up and sort time points
sample_data(phy)$Visit <- sample_data(phy)$Visit %>% 
  recode('Week1'='1w', '1mth'='1m', '1yr'='1y', 
         '4yr'  ='4y', '5yr' ='5y', '6yr'='6y') %>% 
  factor(levels = c('1w', '1m', '1y', '4y','5y', '6y'))

# Only include fecal samples from the child
# Remove rare taxa
phy <- phy %>%
  subset_samples(Sampletype == 'Faeces') %>%
  filter_taxa(function(x) sum(x) > 110, TRUE)




# -------------------------------------------------------------- #
# Characterize sequencing effort and generate rarefaction curves #
# -------------------------------------------------------------- #

# Minimum and mean depths
sample_data(phy)$Reads %>% mean %>% round
sample_data(phy)$Reads %>% min

# Colors for time points
cols <- c('#d94242', '#dc9439', '#249262', '#15638c', '#6b5090', '#c072a0')

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
sum(sample_data(phy)$Reads > 200000)

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
        axis.title.y=element_blank())

# Generate rarefaction curves for all samples ~1h on a good laptop
rarefaction_curves <- ggrare(phy, step = 50, se = FALSE, color = 'Visit')
rarefaction_curves <- rarefaction_curves$data
#load('data/rarefaction_curves.RData')

# Rarefaction curves of all samples
plot_rare_all <- rarefaction_curves %>%
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
plot_rare_separate <- rarefaction_curves %>%
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




# --------------------------------- #
# Alpha diversity using rarefaction #
# --------------------------------- #

# ASV counts in format compatible with vegan
counts <- phy %>% otu_table() %>% as.data.frame() %>% t()

# List of data frames that contain alpha diversity measures for each sample
# Each element is randomly subsampled prior to calculating diversity
alphas <- lapply(1:100, function(i) {
  counts %>%
    vegan::rrarefy(min(rowSums(counts))) %>%
    t() %>% as.data.frame() %>% rownames_to_column(var = "ASV") %>%
    pivot_longer(-ASV, names_to = 'Sample', values_to = 'Reads') %>%
    group_by(Sample) %>%
    summarize(Richness = sum(Reads > 0),
              Shannon = vegan::diversity(Reads, index = 'shannon'),
              Simpson = vegan::diversity(Reads, index = 'simpson'),
              Invsimpson = vegan::diversity(Reads, index = 'invsimpson'),
              Chao1 = estimateR(Reads)[[2]])
  }
)

# Alpha diversity mean across all data frames
alpha <- transpose(alphas)[-1] %>% 
  map(function(x) rowMeans(do.call(rbind.data.frame, transpose(x)))) %>% 
  bind_cols() %>% 
  mutate(.sample = alphas[[1]]$Sample) %>% 
  left_join(as_tibble(sample_data(phy)[, c('ABCNO','Visit','Run','Reads')]), .)
#load('data/alpha.Rdata')

# Simpson diversity distributions by time point
simpson <- alpha %>%
  ggplot(aes(y=Invsimpson, x=Visit)) +
  geom_half_violin(aes(fill=Visit), alpha = 0.4, width=0.6, side='r', scale = "width") +
  geom_half_point(aes(col=Visit), side = 'l', size=0.6, alpha=0.3) +
  geom_boxplot(aes(fill=Visit), outlier.shape=NA, width=0.2, notch=T) +
  theme_bw() +
  ylab('Inverse Simpson Diversity Index') +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  theme(axis.title.x=element_blank(), legend.position='none') +
  stat_compare_means(comparisons = list(c('4y', '5y'), c('4y', '6y'))) +
  stat_compare_means(comparisons = list(c('1w', '1m'))) +
  stat_compare_means(comparisons = list(c('5y', '6y')))

# Shannon diversity distributions by time point
shannon <- alpha %>%
  ggplot(aes(y=Shannon, x=Visit)) +
  geom_half_violin(aes(fill=Visit), alpha = 0.4, width=0.6, side='r', scale = "width") +
  geom_half_point(aes(col=Visit), side = 'l', size=0.6, alpha=0.3) +
  geom_boxplot(aes(fill=Visit), outlier.shape=NA, width=0.2, notch=T) +
  theme_bw() +
  ylab('Shannon Diversity Index') +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  theme(axis.title.x=element_blank(), legend.position='none') +
  stat_compare_means(comparisons = list(c('4y', '5y'), c('4y', '6y'))) +
  stat_compare_means(comparisons = list(c('1w', '1m'))) +
  stat_compare_means(comparisons = list(c('5y', '6y')))

# ASV Richness distributions by time point
richness <- alpha %>%
  ggplot(aes(y=Richness, x=Visit)) +
  geom_half_violin(aes(fill=Visit), alpha = 0.4, width=0.6, side='r', scale = "width") +
  geom_half_point(aes(col=Visit), side = 'l', size=0.6, alpha=0.3) +
  geom_boxplot(aes(fill=Visit), outlier.shape=NA, width=0.2, notch=T) +
  theme_bw() +
  ylab('ASV Richness') +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  theme(axis.title.x=element_blank(), legend.position='none') +
  stat_compare_means(comparisons = list(c('4y', '5y'), c('4y', '6y'))) +
  stat_compare_means(comparisons = list(c('1w', '1m'))) +
  stat_compare_means(comparisons = list(c('5y', '6y')))

# Correlation between Richness and Shannon diversity 
shannon_cor <- alpha %>%
  ggscatter(x='Richness', y='Shannon', col='Visit', fullrange = F, size = 0.85, 
            add='reg.line', alpha=0.1, conf.int=F, palette=cols) +
  stat_cor(aes(color = Visit), method = 'spearman', size = 3,
           label.x.npc = 0.7, label.y.npc = 'bottom') +
  ylab('Shannon Diversity Index') +
  xlab('ASV Richness') +
  theme_bw() +
  theme(legend.position = 'none')

# Correlation between Richness and Simpson diveristy
simpson_cor <- alpha %>%
  ggscatter(x='Richness', y='Invsimpson', col='Visit', fullrange = F, size = 0.85, 
            add='reg.line', alpha=0.1, conf.int=F, palette=cols) +
  stat_cor(aes(color = Visit), method = 'spearman', size = 3,
           label.x.npc = 'left', label.y.npc = 'top') +
  ylab('Inverse Simpson Diversity Index') +
  xlab('ASV Richness') +
  theme_bw() +
  theme(legend.position = 'none')
  
((shannon_cor/simpson_cor) | richness) / (simpson + shannon)
  





# ------------------------------------------------------ #
# Effect of sequencing depth on alpha diversity measures #
# ------------------------------------------------------ #

# Alpha diversity values with and without rarefaction
diversity <- phy %>%
  otu_table() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'ASV') %>% 
  pivot_longer(-ASV, names_to='.sample', values_to='Reads') %>% 
  group_by(.sample) %>% 
  summarize(Invsimpson = diversity(Reads, index='invsimpson'),
            Richness = sum(Reads > 0)) %>% 
  list(Rarefaction = select(alpha, .sample, Richness, Invsimpson), Raw = .) %>%
  imap(~pivot_longer(., -.sample, names_to = 'Diversity', values_to = .y)) %>%
  reduce(full_join) %>% 
  pivot_longer(cols = c(Raw, Rarefaction), names_to = 'Rarefaction', 
               values_to = 'Value') %>% 
  left_join(select(alpha, .sample, ABCNO, Visit, Run, Reads))

# Richness correlation with depth colored by visit
plot_richness <- diversity %>% 
  filter(Diversity == 'Richness') %>% 
  ggscatter(x='Reads', y='Value', col='Visit', size = 0.85,
            add='reg.line', alpha=0.1, conf.int=T, palette=cols) +
  facet_wrap(~Rarefaction) +
  xlim(0, 200000) +
  ylim(0, 400) +
  theme_bw() +
  stat_cor(aes(color = Visit), method = 'spearman', size = 3.5,
           label.x.npc = 'left', label.y.npc = 'top') +
  ylab('ASV Richness') +
  xlab('Sample sequencing depth') +
  theme(legend.position = 'none',
        strip.text = element_text(size = 11),
        strip.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Shannon correlation with depth colored by visit
plot_simpson <- diversity %>% 
  filter(Diversity == 'Invsimpson') %>% 
  ggscatter(x='Reads', y='Value', col='Visit', fullrange = F, size = 0.85, 
            add='reg.line', alpha=0.1, conf.int=T, palette=cols) +
  facet_wrap(~Rarefaction) +
  xlim(0, 200000) +
  ylim(0, 60) +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 50000, 100000, 150000, 200000), 
                     labels = c(0, '50k', '100k', '150k', '200k'),
                     limits = c(0, 200000)) +
  stat_cor(aes(color = Visit), method = 'spearman', size = 3.5,
           label.x.npc = 'left', label.y.npc = 'top') +
  ylab('Inverse Simpson Index') +
  xlab('Sample sequencing depth') +
  theme(legend.position = 'none',
        strip.text.x = element_blank())

# 800x700
plot_richness / plot_simpson





# ------------------------------------------------------- #
# Characterizing depth biases on the run level at 6 years #
# ------------------------------------------------------- #

# Table with information on 6 year runs
Runs <- sample_data(phy) %>% 
  as_tibble() %>% 
  group_by(Run) %>%
  summarize(Samples = n())
Runs <- sample_data(phy) %>% 
  as_tibble() %>% 
  filter(Visit == '6y') %>% 
  group_by(Run) %>%
  summarize(six = n()) %>% 
  left_join(Runs) %>% 
  mutate('Six years' = (six/Samples*100) %>% round() %>% paste0('%')) %>% 
  select(-six)
Runs <- diversity %>% 
  filter(Visit=='6y', Diversity=='Shannon', Rarefaction=='Rarefaction') %>% 
  group_by(Run) %>% 
  cor_test(Reads, Value, method = 'spearman') %>% 
  left_join(Runs, .) %>% 
  select(-var1, -var2, -statistic, -method) %>% 
  rename('Spearman' = 'cor', 'P-value' = 'p')
table <- ggtexttable(Runs, rows = NULL, theme = ttheme(tbody.style = tbody_style(
  color = c('#3A4340', '#f2aed1', '#9A8C9A', '#C92C94'), fill = 0,
  hjust = c(rep(0, nrow(Runs)), rep(1, nrow(Runs)*4)),
  x = c(rep(.1, nrow(Runs)), rep(.9, nrow(Runs)*4))),
  colnames.style = colnames_style(color = 'black', fill = 0, hjust=0, x=0.1)))

# Simpson diversity by depth colored by 6-year run
plot_simpson_6y <- diversity %>% 
  filter(Visit == '6y',
         Diversity == 'Invsimpson',
         Rarefaction == 'Rarefaction') %>% 
  ggscatter(x='Reads', y='Value', col='Run', fullrange = F, 
            add='reg.line', alpha=0.3, conf.int=F) +
  xlim(0, 105000) +
  stat_cor(aes(color = Run), method = 'spearman', size = 3.5,
           label.x.npc = 'left', label.y.npc = 'top') +
  ylab('Inverse Simpson Diversity Index') +
  xlab('Sample Sequencing Depth') +
  theme(legend.position = 'none',
        strip.text.x = element_blank()) +
  scale_color_manual(values = c('#3A4340', '#f2aed1', '#9A8C9A', '#C92C94')) +
  scale_fill_manual(values = c('#3A4340', '#f2aed1', '#9A8C9A', '#C92C94'))



