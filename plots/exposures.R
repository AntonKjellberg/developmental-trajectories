library(tidyverse) # Data processing
library(speedyseq) # Phyloseq with faster functions
library(rstatix) # Pipe-friendly statistic tests
library(patchwork) # Combine plots
library(ggpubr) # Statistics and tables in plots
library(gghalves) # Half-half plots

# Metadata on environmental factors
exposures <- readxl::read_excel('data/exposures.xlsx')
# Copy of the raw data
exposures_copy <- exposures
# Phyloseq object with raw counts and developmental trajectories in sample data
load('data/phy_assigned.RData')
  
# CLR-transform
phy <- phy %>% 
  microbiome::transform('clr')




# -------------------------------------------------------------------- #
# Characterize distributions and confounders in environmental factors  #
# -------------------------------------------------------------------- #

# Custom legend for the next plot
sum(exposures$breast_exclusive < 7, na.rm=T)
legend <- data.frame(x = 1, y = 1, Lines = c('<7 Days (49)','Quartiles')) %>% 
  ggplot(aes(x=x, y=y)) +
  geom_line(aes(color=Lines), linetype='dashed') +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = c('red', 'black')) +
  labs(color = NULL) +
  theme_bw() +
  theme(legend.key = element_rect(fill=NA),
        legend.background = element_rect(fill=NA))

# Breastfeeding distribution
plot_breastfeeding <- exposures %>% 
  ggplot(aes(x=breast_exclusive)) +
  geom_histogram(fill = 'lightgrey', alpha = 0.7, binwidth = 30, boundary = 0, closed = 'left') +
  geom_vline(xintercept = quantile(exposures$breast_exclusive, c(0.25, 0.5, 0.75), na.rm=T), 
             linetype = 'dashed') +
  geom_vline(xintercept = 7, col = 'red', 
             linetype = 'dashed') +
  ylab('Number of children') +
  xlab('Days of exclusive breastfeeding') +
  scale_x_continuous(breaks = c(0, 30, quantile(exposures$breast_exclusive, c(0.25, 0.5, 0.75), na.rm = T), 500, 800), 
                     labels = c(0, 30, quantile(exposures$breast_exclusive, c(0.25, 0.5, 0.75), na.rm = T), 500, 800),
                     limits = c(0,800), minor_breaks = NULL) +
  theme_bw() +
  inset_element(as_ggplot(get_legend(legend)), 1.5,0,0,1.5)

# Distribution of landuse intensity by pets
plot_landuse <- exposures %>% 
  mutate(Pet = ifelse(cat_birth == 1, 'Cat at birth', NA)) %>% filter(!is.na(Pet)) %>% 
  rbind(mutate(exposures, Pet = ifelse(dog_birth == 1, 'Dog at birth', NA)) %>% filter(!is.na(Pet))) %>% 
  rbind(mutate(exposures, Pet = 'All Children')) %>% filter(!is.na(Pet)) %>% 
  filter(landuse_gradient > -100) %>% 
  ggplot(aes(x=Pet, y=landuse_gradient)) +
  geom_half_violin(aes(fill=Pet), alpha = 0.4, side='r', width = 1.1) +
  geom_half_point(aes(col=Pet), side = 'l', size=0.6, alpha=0.3) +
  geom_boxplot(aes(fill=Pet), outlier.shape=NA, width=0.2, notch=T) +
  theme_bw() +
  scale_color_manual(values=c('lightgrey', '#fc926f', '#b5907b')) +
  scale_fill_manual(values=c('lightgrey', '#fc926f', '#b5907b')) +
  theme(axis.title.y=element_blank(), legend.position='none') +
  ylab('Landuse Intensity Index') +
  scale_y_reverse(breaks = c(0, -0.05)) +
  coord_flip()




# -------------------- #
# Make bubblegum plot  #
# -------------------- #

# Clean up tibble with exposures
exposures <- exposures %>% 
  rename(ABCNO = abcno) %>% 
  select(-ABCNO, ABCNO) %>% 
  # mutate(`Breastfeeding 1w` = ifelse(breast_exclusive > 7, 'Yes', 'No')) %>% 
  # mutate(`Breastfeeding 1m` = ifelse(breast_exclusive > 30, 'Yes', 'No')) %>% 
  # mutate(`Breastfeeding 1y` = ifelse(breast_exclusive > 365, 'Yes', 'No')) %>% 
  mutate(`Breastfeeding 6m` = ifelse(breast_exclusive > 180, 'Yes', 'No')) %>%
  mutate(Cesarean = ifelse(delivery %in% c('Acute sectio', 'Planned sectio'), 'Yes', 'No')) %>%
  mutate(Antibiotics = ifelse(ANTIBIOTICSBIRTH != 'No', 'Yes', ANTIBIOTICSBIRTH)) %>%
  mutate(Antibiotics = ifelse(Cesarean == 'Yes', NA, Antibiotics)) %>% 
  mutate(Cesarean = ifelse(Antibiotics == 'Yes' & !is.na(Antibiotics), NA, Cesarean)) %>% 
  mutate(Siblings = ifelse(oldchild == 0, 'No', 'Yes')) %>%
  mutate(Rural = ifelse(landuse_gradient <= quantile(landuse_gradient, 0.25, na.rm=T), 'Yes', 'No')) %>%
  mutate(Rural = ifelse(landuse_gradient >= quantile(landuse_gradient, 0.25, na.rm=T) &
                          landuse_gradient <= quantile(landuse_gradient, 0.75, na.rm=T), NA, Rural)) %>% 
  mutate(Pets = ifelse(catordog == 0, 'No', 'Yes')) %>% 
  mutate(Pets = ifelse(landuse_gradient < quantile(landuse_gradient, 0.75, na.rm=T), NA, Pets)) %>%
  select(ABCNO:length(.))

# Colonizer abundances over time by ABCNO
trajectories <- tax_glom(phy, taxrank = 'Colonizer') %>% 
  psmelt() %>% 
  group_by(ABCNO, Colonizer, Time) %>%
  reframe(Abundance = Abundance) %>% 
  dplyr::left_join(exposures, by = 'ABCNO')

# Colonizer mean abundance differences over time by environmental exposure
bubblegum <- trajectories %>%  
  pivot_longer(cols = tail(names(.), length(exposures)-1), 
               names_to = 'Exposure', 
               values_to = 'Value') %>% 
  mutate(Exposure = as_factor(Exposure)) %>% 
  group_by(Exposure, Value, Time, Colonizer) %>% 
  reframe(Abundance = mean(Abundance), na.omit = T) %>% 
  group_by(Time, Colonizer, Exposure) %>%
  summarize(Difference = diff(Abundance[Value %in% c('Yes', 'No')]))

# Add p-values
bubblegum <- trajectories %>%
  group_by(Time, Colonizer) %>%
  summarise(across(colnames(exposures[-1]), 
                   ~ wilcox.test(Abundance ~ .) %>% tidy())) %>% 
  as.matrix %>% as_tibble() %>% 
  select(1:2, contains('p.value')) %>% 
  rename_with(~gsub('\\.p\\.value', '', .)) %>% 
  pivot_longer(-c(Time, Colonizer), names_to='Exposure', values_to='p') %>%
  left_join(bubblegum) %>% 
  mutate(p = as.numeric(p)) %>% 
  mutate(Label = ifelse(p < 0.05, '*', ''))

# Add number of samples in each group
bubblegum$Exposure <- exposures %>%
  pivot_longer(-1, names_to = 'Exposure', values_to = 'Value') %>%
  mutate(Exposure = as_factor(Exposure)) %>% 
  group_by(Exposure) %>%
  summarize(Yes = sum(Value == 'Yes', na.rm = TRUE),
            No  = sum(Value == 'No', na.rm = TRUE)) %>% 
  mutate(Exposure_percentage = paste0(Exposure, '\n', Yes, ' vs ', No)) %>% 
  .[[4]] %>% rep(5*length(unique(bubblegum$Colonizer))) %>% 
  factor(levels = unique(.))

# Add number of samples at each time point
bubblegum$Time <- trajectories %>% 
  count(Time) %>% 
  mutate(Time_count = paste0(Time, '\n', n/4)) %>% 
  .[[3]] %>% rep(each = (length(exposures)-1)*length(unique(bubblegum$Colonizer))) %>% 
  factor(levels = unique(.))

# Bubblegum plot
# Group p-value size above 0.15 and below 0.001
plot_bubblegum <- bubblegum %>% 
  mutate(p = case_when(p < 0.001 ~ 0.001, p > 0.15 ~ 0.15, TRUE ~ p)) %>% 
  ggplot() + 
  geom_point(aes(x=Time, y=Exposure, size=sqrt(-log(p))-1, fill=Difference), shape=21) +
  facet_wrap(~Colonizer, scale = 'free_y') +
  scale_fill_gradient2(low = 'darkblue', mid = 'white', high = 'darkred', 
                       name='Change') +
  geom_text(mapping=aes(x=Time, y=Exposure, label=Label), size=3) + 
  scale_y_discrete(limits=rev) +
  scale_size_area(max_size = 11, name='P-value',
                  breaks = c(sqrt(-log(0.001))-1,
                             sqrt(-log(0.01))-1,
                             sqrt(-log(0.05))-1,
                             sqrt(-log(0.15))-1),
                  labels = c('< 0.001', '   0.01', '   0.05', '> 0.15')) +
  theme_bw() +
  theme(strip.background = element_rect(colour=NA, fill=NA),
        aspect.ratio = 0.7,
        legend.position='bottom', 
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

# 6.25x6.25
plot_bubblegum
# 15x9
(plot_landuse + plot_breastfeeding) / plot_bubblegum +
  plot_layout(heights = c(1, 4))




# --------------------------------------------------------------- #
# Characterize Rumminococcus gnavus in the context of antibiotics #
# --------------------------------------------------------------- #

# Alpha diversity data in each sample
load('data/alpha.Rdata')
# Colors for time points
cols <- c('#d94242', '#dc9439', '#249262', '#15638c', '#6b5090', '#c072a0')

# Ruminococcus gnavus abundance by antibiotics at birth at each time point
gnavus1 <- phy %>%
  psmelt() %>%
  filter(Species == 'Ruminococcus_B gnavus') %>%
  select(ABCNO, Abundance, Time) %>%
  left_join(exposures_copy %>% rename(ABCNO = abcno)) %>%
  mutate(Antibiotics = ifelse(ANTIBIOTICSBIRTH != 'No', 'Yes', ANTIBIOTICSBIRTH)) %>%
  ggplot(aes(y=Abundance, x=Antibiotics)) +
  facet_wrap(~Time, nrow = 1) +
  geom_half_violin(aes(fill=Time), alpha = 0.4, side='r', width = 1, scale = "width") +
  geom_half_point(aes(col=Time), side = 'l', size=0.6, alpha=0.3) +
  geom_boxplot(aes(fill=Time), outlier.shape=NA, width=0.15) +
  stat_compare_means(comparisons = list(c('Yes', 'No'))) +
  scale_fill_manual(values=cols[-5]) +
  scale_color_manual(values=cols[-5]) +
  theme_bw() +
  theme(legend.position='none',
        strip.background = element_rect(colour=NA, fill=NA)) +
  xlab('Antibiotics at birth') +
  ylab('Rumminococcus Gnavus CLR-Transformed Abundance') 

# Shannon diveristy by antibiotics at birth at each time point
gnavus2 <- phy %>% 
  psmelt() %>%
  filter(Species == 'Ruminococcus_B gnavus') %>% 
  left_join(alpha, by = c('ABCNO', 'Visit')) %>% 
  left_join(exposures_copy %>% rename(ABCNO = abcno), by = 'ABCNO') %>% 
  mutate(Antibiotics = ifelse(ANTIBIOTICSBIRTH != 'No', 'Yes', ANTIBIOTICSBIRTH)) %>% 
  ggplot(aes(y=Shannon, x=Antibiotics)) +
  facet_wrap(~Time, nrow=1) +
  geom_half_violin(aes(fill=Time), alpha = 0.4, side='r', width = 1, scale = "width") +
  geom_half_point(aes(col=Time), side = 'l', size=0.6, alpha=0.3) +
  geom_boxplot(aes(fill=Time), outlier.shape=NA, width=0.15) +
  stat_compare_means(comparisons = list(c('Yes', 'No'))) +
  theme_bw() +
  theme(legend.position='none',
        strip.background = element_rect(colour=NA, fill=NA)) +
  scale_fill_manual(values=cols[-5]) +
  scale_color_manual(values=cols[-5]) +
  xlab('Antibiotics at birth') +
  ylab('Shannon Diversity Index') 

# Shannon index by Ruminococcus gnavus abundance
gnavus3 <- phy %>% 
  psmelt() %>%
  filter(Species == 'Ruminococcus_B gnavus') %>% 
  left_join(alpha, by=c('ABCNO', 'Visit')) %>% 
  ggscatter(x='Abundance', y='Shannon', col='Visit', fullrange=F, size=1, 
            add='reg.line', alpha=0.17, conf.int=T, palette=cols) +
  theme_bw() +
  stat_cor(aes(color = Visit), method = 'pearson', size = 3.5,
           label.x.npc = 0.765, label.y.npc = 1) +
  ylab('Shannon Index') +
  xlab('Ruminococcus gnavus CLR-Transformed Abundance') +
  theme(legend.position = 'none',
        strip.text.x = element_blank())

# 13x7
gnavus1 / gnavus2 / gnavus3 + plot_layout(heights = c(1, 1, 0.6))

