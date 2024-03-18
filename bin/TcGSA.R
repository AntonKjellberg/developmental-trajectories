library(phyloseq)
library(tidyverse)
library(microbiome)
library(TcGSA)
library(patchwork)

setwd('~/R/copsac/github')
load('data/dada2_gtdb_2021_11_18.RData')
load('data/Feaces_16S_Traj_DMM.RData')

# Split 5y into 4y and 6y. 4y if no sample at 4y already. Otherwise 5y
sample_data(phy) <- phy %>% 
  sample_data %>% 
  data.frame %>% 
  mutate(Time = case_when(
    Time == "5y" & ABCNO %in% filter(., Time == '4y')$ABCNO ~ "6y",
    Time == "5y" ~ "4y",
    TRUE ~ Time))
# Clean up time points
sample_data(phy)$Visit <- sample_data(phy)$Visit %>% 
  recode('Week1'='1w', '1mth'='1m', '1yr'='1y', 
         '4yr'  ='4y', '5yr' ='5y', '6yr'='6y') %>% 
  factor(levels = c('1w', '1m', '1y', '4y', '5y', '6y'))
sample_data(phy)$Time <- sample_data(phy)$Time %>% 
  factor(levels = c('1w', '1m', '1y', '4y', '6y'))

# Only include fecal samples from the child
# Remove taxa that are not present
phy <- phy %>% 
  subset_samples(Visit %in% c('1w', '1m', '1y', '4y', '5y', '6y')) %>%
  subset_samples(Sampletype == 'Faeces')

# Based on relative abundance
phy <- speedyseq::tax_glom(phy, "Genus") %>% microbiome::transform('compositional')

# Subset of genera present in more than 20% of samples
# With a mean relative abundance above 0.001% of samples at a time point
# The 5 year time point was not split up for this but was for clustering
genera <- c()
for (visit in levels(sample_data(phy)$Visit)) {
  subset <- phy %>% 
    subset_samples(Visit == visit) %>% 
    filter_taxa(function(x) sum(x > 0) / length(x) >= 0.2, TRUE)
  # Calculate mean relative abundance
  subset <- data.frame(Genus = as.character(tax_table(subset)[,"Genus"]),
                       mra = taxa_sums(subset) / nsamples(subset))
  genera <- c(genera, filter(subset, mra > 1e-3) %>% pull(Genus)) %>% unique
}

# Filter out uninteresting genera
phy <- prune_taxa(tax_table(phy)[,"Genus"] %in% genera, phy)

# Extract data for the TcGSA clustering
# Genera mean relative abundances were log transformed
# A pseudo count corresponding to the mean relative abundance cutoff was added
#
# I've tried combinations of Z-score, CLR, Log10, MRA, ranks, ect.
# microbiome::transform('log10p')
# microbiome::transform('Z')
# microbiome::transform('clr')
otu_table <- phy %>% 
  transform_sample_counts(function(x) log(1e-3+x)) %>% 
  otu_table() %>% 
  data.frame(row.names = as_tibble(tax_table(phy))$Genus)

genera <- list(genesets = list(as_tibble(tax_table(phy))$Genus),
               geneset.names = "AllGenera",
               geneset.descriptions = "All genera from the dataset")
IDs <- sample_data(phy) %>% row.names()
# I've tried scaling the time points differently
# 1w and 1m time points are swapped in the original script
# I decided to correct it here. No influence on the clustering
time <- sample_data(phy)$Visit %>% recode('1w'=1,'1m'=2,'1y'=3,
                                          '4y'=4,'5y'=5, '6y'=6)

# Use TcGSA to assign colonizers
trajectories <- plot1GS(expr=otu_table, gmt=genera,
                        Subject_ID = IDs, TimePoint = time,
                        geneset.name="AllGenera",
                        aggreg.fun = "mean")

# Plot the estimates and means
trajectories$p$data %>%
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
               fun.y="median", geom="line", col="black", linewidth=0.9) +
  theme_minimal() +
  scale_color_manual(values=c('darkgreen','orange', 'darkred', 'darkblue'))

# Add colonizer assingments to the phyloseq object
tax_table(phy) <- tax_table(phy) %>% 
  data.frame() %>% 
  dplyr::left_join(trajectories$classif, by = c("Genus" = "ProbeID")) %>%
  data.frame(row.names = rownames(tax_table(phy))) %>% 
  dplyr::rename(Colonizer = 'Cluster') %>% 
  relocate(Colonizer, .before=Domain) %>% 
  as.matrix() %>%
  tax_table

# Plot colonizer means
tax_glom(phy, taxrank = "Colonizer") %>%
  speedyseq::psmelt() %>%
  group_by(Colonizer, Visit) %>%
  summarize(Mean = mean(Abundance, na.rm = TRUE)) %>%
  ggplot(aes(y=Mean, x=Visit, group=Colonizer, col=Colonizer)) +
  geom_line(linewidth=1) +
  theme_bw() +
  scale_color_manual(values=c('darkgreen','orange', 'darkred', 'darkblue')) +
  ylab('Mean Relative Abundance')

# Plot 1
plot1 <- subset_taxa(phy, Colonizer == 1) %>%
  filter_taxa(function(x) mean(x) > 0.00145, T) %>% 
  psmelt() %>%
  group_by(Genus, Visit) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ggplot(aes(y=Abundance, x=Visit, group=Genus, col=Genus, linetype=Genus)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  ylab('Mean Relative Abundance') +
  scale_color_manual(name = 'Early ASVs', 
                     values = rep(c('darkgreen', '#89ab41', '#72C79B'), 4)) +
  scale_linetype_manual(name = 'Early ASVs', 
                        values = rep(c('solid','dashed','dotted','dotdash'), each=3)) +
  theme(legend.title = element_text(face = 'bold'),
        axis.title.x=element_blank())


# Plot 4
plot4 <- subset_taxa(phy, Colonizer == 4) %>%
  filter_taxa(function(x) mean(x) > 0.00145, T) %>% 
  psmelt() %>%
  group_by(Genus, Visit) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ggplot(aes(y=Abundance, x=Visit, group=Genus, col=Genus, linetype=Genus)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  ylab('Mean Relative Abundance') +
  scale_color_manual(name = 'Cluster 4 ASVs', 
                     values = rep(c('orange', '#B69006', '#FFDB58'), 4)) +
  scale_linetype_manual(name = 'Cluster 4 ASVs', 
                        values = rep(c('solid','dashed','dotted','dotdash'), each=3)) +
  theme(legend.title = element_text(face = 'bold'),
        axis.title.x=element_blank())

# Plot 3
plot3 <- subset_taxa(phy, Colonizer == 3) %>%
  filter_taxa(function(x) mean(x) > 0.0001, T) %>% 
  psmelt() %>%
  group_by(Genus, Visit) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ggplot(aes(y=Abundance, x=Visit, group=Genus, col=Genus, linetype=Genus)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  ylab('Mean Relative Abundance') +
  scale_color_manual(name = 'Cluster 3 ASVs', 
                     values = rep(c('darkred', 'red', '#d16680'), 4)) +
  scale_linetype_manual(name = 'Cluster 3 ASVs', 
                        values = rep(c('solid','dashed','dotted','dotdash'), each=3)) +
  theme(legend.title = element_text(face = 'bold'),
        axis.title.x=element_blank())

# Plot 2
plot2 <- subset_taxa(phy, Colonizer == 2) %>%
  filter_taxa(function(x) mean(x) > 0.0058, T) %>% 
  psmelt() %>%
  group_by(Genus, Visit) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ggplot(aes(y=Abundance, x=Visit, group=Genus, col=Genus, linetype=Genus)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  ylab('Mean Relative Abundance') +
  scale_color_manual(name = 'Cluster 2 ASVs', 
                     values = rep(c('darkblue', 'cyan4', 'lightblue'), 4)) +
  scale_linetype_manual(name = 'Cluster 2 ASVs', 
                        values = rep(c('solid','dashed','dotted','dotdash'), each=3)) +
  theme(legend.title = element_text(face = 'bold'),
        axis.title.x=element_blank())

# 700x1400
plot1 + plot4 + plot3 + plot2 + 
  plot_layout(ncol = 1) & theme(legend.justification = 'left')
