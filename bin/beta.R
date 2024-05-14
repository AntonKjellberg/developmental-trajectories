library(tidyverse) # Data processing
library(speedyseq) # Phyloseq with faster functions
library(vegan) # Community ecology processing
library(usedist) # Distance matrix utilities
library(ggside) # Plots on the sides of plots
library(ggpubr) # Plotting and statistics

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
# Remove bleed-through taxa
phy <- phy %>% 
  subset_samples(Visit %in% c('1w', '1m', '1y', '4y', '5y', '6y')) %>%
  subset_samples(Sampletype == 'Faeces') %>% 
  filter_taxa(function(x) sum(x) > 110, TRUE)




# ------------------------------------------------------------ #
# Characterize beta diversity using rarefaction and NMDS plots #
# ------------------------------------------------------------ #

# Matrix with otu counts
counts <- phy %>%
  otu_table() %>% 
  as.data.frame() %>% 
  t()

# Create a Bray Curtis dissimilarity matrix with rarefaction
# 100 iterations of subsampling to 2016 read takes 10 hours on a modern laptop
dist <- avgdist(counts, sample = min(rowSums(counts)), iterations = 100)
#load('data/dist_rare_2016_100.Rdata')

# Distance ordination using Non-metric Multidimensional Scaling
nmds <- metaMDS(dist)

# Colors for the time points
cols <- c('#d94242', '#dc9439', '#249262', '#15638c', '#6b5090', '#c072a0')

# NMDS plot with all time points
plot_all <- sample_data(phy) %>% 
  as_tibble() %>% 
  bind_cols(as_tibble(scores(nmds))) %>% 
  ggscatter(x='NMDS1', y='NMDS2', col='Visit', alpha=0.4, ellipse = T, 
            ellipse.alpha = 0) +
  xlim(c(-0.52, 0.52)) +
  ylim(c(-0.52, 0.52)) +
  geom_ysideboxplot(aes(x = Visit, fill = Visit, col = Visit), alpha = 0.4,
                    orientation = "x",  outlier.shape = NA, notch = T) +
  geom_xsideboxplot(aes(y = Visit, fill = Visit, col = Visit), alpha = 0.4,
                    orientation = "y",  outlier.shape = NA, notch = T) +
  theme_ggside_void() +
  theme(ggside.axis.line.x = element_blank(), 
        ggside.panel.scale = 0.2,
        legend.position = 'none',
        axis.title.x = element_text(hjust = 0.41),
        axis.title.y = element_text(hjust = 0.41)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols)

# Subset the distance matrix to only include late samples and perform NMDS
# Subsetting the NMDS scores directly is not optimal 
nmds_late <- sample_data(phy) %>% 
  as_tibble() %>% 
  filter(Visit %in% c('4y', '5y', '6y')) %>% 
  .$SampleID %>% 
  dist_subset(dist, .) %>% 
  metaMDS()

# NMDS plot with late samples
plot_late <- sample_data(phy) %>% 
  as_tibble() %>% 
  filter(Visit %in% c('4y', '5y', '6y')) %>% 
  bind_cols(as_tibble(scores(nmds_late))) %>% 
  ggscatter(x='NMDS1', y='NMDS2', col='Visit', alpha=0.4, ellipse = T, 
            ellipse.alpha = 0) +
  xlim(c(-0.52, 0.52)) +
  ylim(c(-0.52, 0.52)) +
  geom_ysideboxplot(aes(x = Visit, fill = Visit, col = Visit), alpha = 0.4,
                    orientation = "x",  outlier.shape = NA, notch = T) +
  geom_xsideboxplot(aes(y = Visit, fill = Visit, col = Visit), alpha = 0.4,
                    orientation = "y",  outlier.shape = NA, notch = T) +
  theme_ggside_void() +
  theme(ggside.axis.line.x = element_blank(), 
        ggside.panel.scale = 0.2,
        legend.position = 'none',
        axis.title.x = element_text(hjust = 0.41),
        axis.title.y = element_text(hjust = 0.41)) +
  scale_color_manual(values = cols[4:6]) +
  scale_fill_manual(values = cols[4:6])

# Subset the distance matrix to include early samples and perform NMDS
# Subsetting the NMDS scores directly is not optimal 
nmds_early <- sample_data(phy) %>% 
  as_tibble() %>% 
  filter(Visit %in% c('1w', '1m')) %>% 
  .$SampleID %>% 
  dist_subset(dist, .) %>% 
  metaMDS()

# NMDS plot with early samples
plot_early <- sample_data(phy) %>% 
  as_tibble() %>% 
  filter(Visit %in% c('1w', '1m')) %>% 
  bind_cols(as_tibble(scores(nmds_early))) %>% 
  ggscatter(x='NMDS1', y='NMDS2', col='Visit', alpha=0.4, ellipse = T, 
            ellipse.alpha = 0) +
  xlim(c(-0.52, 0.52)) +
  ylim(c(-0.52, 0.52)) +
  geom_ysideboxplot(aes(x = Visit, fill = Visit, col = Visit), alpha = 0.4,
                    orientation = "x",  outlier.shape = NA, notch = T) +
  geom_xsideboxplot(aes(y = Visit, fill = Visit, col = Visit), alpha = 0.4,
                    orientation = "y",  outlier.shape = NA, notch = T) +
  theme_ggside_void() +
  theme(ggside.axis.line.x = element_blank(), 
        ggside.panel.scale = 0.2,
        legend.position = 'none',
        axis.title.x = element_text(hjust = 0.41),
        axis.title.y = element_text(hjust = 0.41)) +
  scale_color_manual(values = cols[1:2]) +
  scale_fill_manual(values = cols[1:2])

plot_all / plot_early / plot_late




# -------------------------------------- #
# Characterize bias in the six-year runs #
# -------------------------------------- #

# Subset the distance matrix to include 6y samples and perform NMDS
nmds_6y <- sample_data(phy) %>% 
  as_tibble() %>% 
  filter(Visit == '6y') %>% 
  .$SampleID %>%
  dist_subset(dist, .) %>% 
  metaMDS()

# NMDS plot with 6y samples colored by sequencing run
plot_6y <- sample_data(phy) %>% 
  as_tibble() %>% 
  filter(Visit == '6y') %>% 
  bind_cols(as_tibble(scores(nmds_6y))) %>% 
  ggscatter(x='NMDS1', y='NMDS2', col='Run', alpha=0.4, ellipse = T, 
            ellipse.alpha = 0) +
  xlim(c(-0.52, 0.52)) +
  ylim(c(-0.52, 0.52)) +
  geom_ysideboxplot(aes(x = Run, fill = Run, col = Run), alpha = 0.4,
                    orientation = "x",  outlier.shape = NA, notch = T) +
  geom_xsideboxplot(aes(y = Run, fill = Run, col = Run), alpha = 0.4,
                    orientation = "y",  outlier.shape = NA, notch = T) +
  theme_ggside_void() +
  theme(legend.position = 'none',
        ggside.axis.line.x = element_blank(), 
        ggside.panel.scale = 0.2,
        axis.title.x = element_text(hjust = 0.41),
        axis.title.y = element_text(hjust = 0.41)) +
  scale_color_manual(values = c('#3A4340', '#f2aed1', '#9A8C9A', '#C92C94')) +
  scale_fill_manual(values = c('#3A4340', '#f2aed1', '#9A8C9A', '#C92C94'))




# ------------------ #
# Pairwise PERMANOVA #
# ------------------ #

# Function to perform pairwise adonis
# phy is the phyloseq object
# group argument is the sample_data column in which all nominal categories are compared
# dist is the distance matrix
# perm is the number of permutations
adonis_pairwise <- function(group, phy, dist, perm=10000) {

# Generate list of pairs to compare
comparisons <- sample_data(phy) %>%
    as_tibble() %>% 
    pull(!!enquo(group)) %>% 
    unique() %>%
    as.character() %>% 
    combn(2, simplify=F)

# perform PERMANOVA on each pair
results <- lapply(comparisons, function(x)
    sample_data(phy) %>% 
    as_tibble() %>% 
    rename(group = !!enquo(group)) %>% 
    filter(group %in% x) %>% 
    adonis2(dist_subset(dist, .$SampleID) ~ group, data = ., permutations = perm) %>% 
    pull('Pr(>F)', 'F') %>% .[1])

# clean up results and fdr adjust
tibble(comparison = comparisons,`F` = results %>% unlist %>% names %>% as.numeric,
       p = results %>% as.numeric) %>%
  mutate(comparison = gsub("[(),\"]", "", comparison)) %>% 
  separate(comparison, into = c("comp1", "comp2"), sep = " ") %>% 
  mutate(comp1 = str_sub(comp1, start = 2)) %>% 
  mutate(p_adj = p.adjust(p, method = 'BH'))

}

# test for differences between time-points
adonis_pairwise(Visit, phy = phy, dist = dist)

# test for differences between six-year runs
adonis_pairwise(Run, phy = subset_samples(phy, Visit == '6y'), dist = dist)



