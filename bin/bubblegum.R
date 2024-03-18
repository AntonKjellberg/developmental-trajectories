library('tidyverse')
library('phyloseq') # I use speedyseq for slow phyloseq functions

setwd('~/R/copsac/github')

load('data/phy_assigned.RData')
load('data/exposure5.RData')
load('data/bmi_6y.RData')

# Clean up the exposures data frame. Make it Yes/No
exposures <- exposures5 %>%
  rename(ABCNO = abcno) %>% 
  dplyr::left_join(bmi_6y[, c(1, 40, 55)], by = 'ABCNO') %>% 
  mutate(Cesarean = ifelse(Delivery %in% c('Acute sectio', 'Planned sectio'), 'Yes', 'No')) %>%
  mutate(Antibiotics = ifelse(ANTIBIOTICSBIRTH != 'No', 'Yes', ANTIBIOTICSBIRTH)) %>%
  mutate(Siblings = ifelse(oldchild == 0, 'No', 'Yes')) %>%
  mutate(Rural = ifelse(urbanrural == 'Urban', 'No', 'Yes')) %>%
  .[, c(1, 9:12)]

# Colonizer abundances by ABCNO per Visit point
trajectories <- tax_glom(phy, taxrank = "Colonizer") %>% 
  speedyseq::psmelt() %>% 
  group_by(ABCNO, Colonizer, Visit) %>%
  reframe(Abundance = Abundance) %>% 
  dplyr::left_join(exposures, by = 'ABCNO')

# Bubblegum data
# Mean abundance differences and minimum group sizes by exposure per Visit point
bubblegum <- tax_glom(phy, taxrank = "Colonizer") %>% 
  speedyseq::psmelt() %>% 
  group_by(ABCNO, Colonizer, Visit) %>%
  reframe(Abundance = Abundance) %>%
  dplyr::left_join(exposures, by = 'ABCNO') %>% 
  pivot_longer(cols = tail(names(.), length(exposures)-1), 
               names_to = "Exposure", 
               values_to = "Value") %>% 
  group_by(Exposure, Value, Visit, Colonizer) %>% 
  reframe(Abundance = mean(Abundance), na.omit = T) %>% 
  group_by(Visit, Exposure, Colonizer) %>%
  summarize(Difference = diff(Abundance[Value %in% c("Yes", "No")]))

# Generate p-values for each bubble using wilcoxon ranked sum test
p_value <- c()
for (visit in unique(bubblegum$Visit)) {
  for (exposure in unique(bubblegum$Exposure)) {
    for (trajectory in unique(bubblegum$Colonizer)) {
      p <- wilcox.test(filter(trajectories, Colonizer == trajectory,
                              !!sym(exposure) == 'No',
                              Visit == visit)$Abundance,
                       filter(trajectories, Colonizer == trajectory,
                              !!sym(exposure) == 'Yes',
                              Visit == visit)$Abundance
      )$p.value
      p_value <- c(p_value, p)
    }
  }
}

# Add p-values to the bubblegum data and label significance
bubblegum$p_value <- p_value
bubblegum <- bubblegum %>% mutate(label = ifelse(p_value < 0.05, "*", ""))

# Add the percentage of children exposed
bubblegum$Exposure <- exposures %>%
  pivot_longer(-1, names_to = 'Exposure', values_to = 'Value') %>%
  group_by(Exposure) %>%
  summarize(Percentage = round(mean(Value == 'Yes', na.rm = TRUE) * 100)) %>%
  mutate(Exposure_percentage = paste0(Exposure, '\n', Percentage, "%")) %>%
  .[[3]] %>% rep(each = 4) %>% rep(6)

# Add the number of samples from each visit
bubblegum$Visit <- trajectories %>% 
  count(Visit) %>% 
  mutate(Visit_count = paste0(Visit, '\n', n/4)) %>% 
  .[[3]] %>% rep(each = 16)

# Bubblegum plot
# Group values above 0.15 and below 0.001
bubblegum %>% 
  mutate(p_value = case_when(
    p_value < 0.001 ~ 0.001,
    p_value > 0.15 ~ 0.15,
    TRUE ~ p_value)) %>% 
  ggplot() + 
  geom_point(aes(x=Visit, y=Exposure, size=sqrt(-log(p_value))-1, fill=Difference), shape=21) +
  facet_wrap(~Colonizer, scale = 'free_y') +
  scale_fill_gradient2(low = 'darkblue', mid = 'white', high = 'darkred', 
                       name='Change') +
  geom_text(mapping=aes(x=Visit, y=Exposure, label=label), size=3) + 
  theme_minimal() + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_size_area(max_size = 11, name='P-value',
                  breaks = c(sqrt(-log(0.001))-1,
                             sqrt(-log(0.01))-1,
                             sqrt(-log(0.05))-1,
                             sqrt(-log(0.15))-1),
                  labels = c('< 0.001', '   0.01', '   0.05', '> 0.15'))

