library(tidyverse) # Data processing
library(speedyseq) # Phyloseq with faster functions
library(ggpubr) # Statistics and tables in plots
library(gghalves) # Half-half plots

# Metadata on environmental factors
exposures <- readxl::read_excel('data/exposures.xlsx')
# Phyloseq object with raw counts and developmental trajectories in sample data
load('data/phy_assigned.RData')

# CLR-transform
phy <- phy %>% 
  microbiome::transform('clr')


# Generate scores
# Only include children with samples from "1w", "1m", "1y" and one of "4y", "5y", "6y"
# Mean of deviation, abundance, and overall abundance mean is taken across "4y", "5y", "6y" if multiple are present
# Fluctuation denotes difference between the abundance change between two time points and the change in mean abundance between the same time points
trajectories <- tax_glom(phy, taxrank = 'Colonizer') %>%
  psmelt() %>% 
  select(ABCNO, Visit, Colonizer, Abundance) %>% 
  group_by(ABCNO) %>%
  filter(any(Visit %in% c("1w", "1m", "1y")) & any(Visit %in% c("4y", "5y", "6y"))) %>% 
  group_by(Colonizer, Visit) %>%
  mutate(Mean = mean(Abundance, na.rm = TRUE)) %>% 
  mutate(Deviation = Abundance - Mean) %>% 
  mutate(Visit = fct_recode(Visit, "4y" = "6y")) %>% 
  mutate(Visit = fct_recode(Visit, "4y" = "5y")) %>% 
  group_by(Colonizer, ABCNO) %>%
  mutate(Deviation = ifelse(Visit == "4y", mean(Deviation[Visit == "4y"]), Deviation)) %>%
  mutate(Mean = ifelse(Visit == "4y", mean(Mean[Visit == "4y"]), Mean)) %>% 
  mutate(Abundance = ifelse(Visit == "4y", mean(Abundance[Visit == "4y"]), Abundance)) %>% 
  rename(Time = 'Visit') %>% 
  group_by(Colonizer, ABCNO) %>%
  distinct %>% 
  ungroup() %>% 
  arrange(ABCNO,  Colonizer, Time) %>% 
  mutate(Fluctuation = c(NA, diff(Abundance - Mean))) %>% 
  mutate(Fluctuation = ifelse(Time == '1w', NA, Fluctuation)) %>% 
  pivot_longer(c('Deviation','Fluctuation'), names_to='Scoring', values_to='Score')

# Plot to see if there is a difference in the mean score across all time points
trajectories %>% 
  group_by(ABCNO, Scoring) %>% 
  summarize(Score = mean(abs(Score), na.rm=T)) %>% 
  left_join(exposures %>% rename(ABCNO = abcno)) %>% 
  mutate(Antibiotics = ifelse(ANTIBIOTICSBIRTH != 'No', 'Yes', ANTIBIOTICSBIRTH)) %>% 
  ggplot(aes(y=Score, x = Antibiotics)) +
  facet_wrap(~Scoring) +
  geom_half_violin(aes(fill=Antibiotics), alpha = 0.4, side='r', width = 1.1, scale = "width") +
  geom_half_point(aes(col=Antibiotics), side = 'l', size=0.6, alpha=0.3) +
  geom_boxplot(aes(fill=Antibiotics), outlier.shape=NA, width=0.15) +
  stat_compare_means(comparisons = list(c('Yes', 'No'))) +
  scale_color_manual(values=c('grey50', 'grey30')) +
  scale_fill_manual(values=c('grey50', 'grey30')) +
  theme_bw() +
  theme(legend.position='none') +
  xlab('Antibiotics at birth') +
  ylab('Score')


