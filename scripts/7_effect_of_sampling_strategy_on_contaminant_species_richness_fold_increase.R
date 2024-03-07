### 7_effect_of_sampling_strategy_on_contaminant_species_richness
### 3. fold increase from one to two samples per pot

# we want to investigate the effect of sampling effort on the observed species 
# richness of contaminant plants and fungi identified by metabarcoding. 

# for this last analysis, we will look into the value of taking a second sample
# per pot. specifically, we will calculate the fold increase in species richness
# between two and one sample. 

### PLANTS
twice <- subset_samples(plants, interval == "double"); twice # subset pots sampled twice
twice <- prune_taxa(taxa_sums(twice)>0, twice); twice # prune taxa with 0 abundance
sort(sample_sums(twice))
twice <- rarefy_even_depth(twice, 
                           replace = F, 
                           trimOTUs = T, 
                           rngseed = 123, 
                           sample.size = 8629); twice

# p46 lost one sample. remove the other one: 
twice <- subset_samples(twice, pot_id!="p46"); twice
#
first <- subset_samples(twice, grepl(pattern = "_1", 
                                    x = twice@sam_data$sample_id)); first = prune_taxa(taxa_sums(first)>0,
                                                                                       first); first
first@sam_data$turn <- "first"; head(first@sam_data) # add a variable to distinguish between first and merged samples
second <- merge_samples(twice, group = "pot_id"); second 

# reassign sample data in second:
second@sam_data$pot_id <- rownames(second@sam_data) # fix pot_id variable
second@sam_data$pot_id <- sub("p","", second@sam_data$pot_id) 
second@sam_data$pot_id <- as.numeric(second@sam_data$pot_id) # this will let us sort by pot_id
second.sam <- data.frame(second@sam_data) %>%  # using arrange() in dplyr
  arrange(pot_id) %>%
  mutate(pot_id=paste0("p",pot_id)); head(second.sam); tail(second.sam)

# now make an identical dataframe to get variables from: 
sam <- data.frame(twice@sam_data) 
sam <- sam %>%
  filter(grepl(pattern = "_2", 
               sam$sample_id))

# and reassign variables: 
second.sam$host <- sam$host
second.sam$container <- sam$container
second.sam$interval <- sam$interval
second.sam$turn <- "second"; # add a variable to distinguish between merged and first samples

# convert til sample_data() object:
second.sam <- sample_data(second.sam); sample_names(second.sam)

# and reassign to twice.m phyloseq object: 
sample_data(second) <- second.sam; twice.m; head(second@sam_data)

# now we can merge the two phyloseq objects in order to normalize them:
fiphy <- merge_phyloseq(first, second); fiphy

# now we have to separate the new object into two dataframes before we combine
# them to calculate and plot fold increase: 
first <- subset_samples(fiphy, turn=="first")
first.df <- data.frame(cbind(estimate_richness(first, 
                                              measures = "Observed"),
                            first@sam_data)) %>% 
  rename("obs_first"="Observed",
         "pot_id_first"="pot_id",
         "turn_first"="turn") %>% 
  mutate(pot_id=pot_id_first); head(first.df)

second <- subset_samples(fiphy, turn=="second")
second.df <- data.frame(cbind(estimate_richness(second, 
                                               measures = "Observed"),
                             second@sam_data)) %>% 
  rename("obs_second"="Observed",
         "pot_id_second"="pot_id",
         "turn_second"="turn") %>% 
  mutate(pot_id=pot_id_second); head(second.df)

fidf <- full_join(first.df, second.df, by = "pot_id") %>% 
  select("obs_first", "obs_second", 
         "pot_id", "pot_id_first", "pot_id_second", 
         "turn_first", "turn_second")

fidf <- fidf %>% 
  mutate(foldincrease = obs_second/obs_first)
ggplot(fidf, aes(x=reorder(pot_id, -foldincrease), 
                          y=foldincrease,
                          fill="chartreuse1", color="black")) + 
  geom_hline(yintercept = c(1, 1.25, 1.5, 1.75, 
                            2, 2.25, 2.5), 
             linetype = 2, color = "grey50") +
  geom_col(width=0.8, linewidth = 0.75) +
  scale_fill_manual(values = "chartreuse1") +
  scale_color_manual(values = "black") +
  theme_classic() +
  theme(legend.position = "none", 
        axis.ticks.x = element_blank()) + 
  ylab("Fold increase") +
  xlab("Sampled pots") +
  scale_x_discrete(labels = paste0(rep(" ", times=20))) +
  labs(title = "A") +
  scale_y_continuous(breaks = seq(0, 3, 0.25), limits = c(0, 3)) +
  coord_cartesian(ylim = c(1, 2.1)) -> fi.p; fi.p

# additionally, we will investigate whether there is a correlation between 
# the richness in the second sample and fold increase for each pot:
cormod <- cor.test(fidf$foldincrease, fidf$obs_second, 
                  method = "pearson"); cormod
ggplot(fidf, aes(x=obs_second, 
                 y=foldincrease,
                 fill="chartreuse1")) + 
  geom_hline(yintercept = c(1, 1.25, 1.5, 1.75, 
                            2, 2.25, 2.5), 
             linetype = 2, color = "grey50") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values = "black") + 
  scale_fill_manual(values = "chartreuse1") +
  geom_point(size=4, 
             stroke = 1.25,
             pch=21, 
             alpha = 1, 
             #position=position_jitter(h=0.01,w=0.01)
             ) + 
  geom_smooth(method = "lm", 
              color="black", 
              fill="grey90", 
              linetype = 1) +
  scale_x_continuous(breaks = seq(0, 20, 5)) +
  xlab("Richness in both samples") +
  ylab(" ") + 
  labs(title = "B") +
  scale_y_continuous(breaks = seq(0, 2.5, 0.25), limits = c(0, 2.5)) +
  coord_cartesian(ylim = c(1, 2.1)) -> fi.point.p; fi.point.p

### FUNGI

twice <- subset_samples(fungi, interval == "double"); twice # subset pots sampled twice
twice <- prune_taxa(taxa_sums(twice)>0, twice); twice # prune taxa with 0 abundance
sort(sample_sums(twice))
twice <- rarefy_even_depth(twice, 
                           replace = F, 
                           trimOTUs = T, 
                           rngseed = 123, 
                           sample.size = min(sample_sums(twice))); twice

first <- subset_samples(twice, grepl(pattern = "_1", 
                                     x = twice@sam_data$sample_id)); first = prune_taxa(taxa_sums(first)>0,
                                                                                        first); first
first@sam_data$turn <- "first"; head(first@sam_data) # add a variable to distinguish between first and merged samples
second <- merge_samples(twice, group = "pot_id"); second 

# reassign sample data in second:
second@sam_data$pot_id <- rownames(second@sam_data) # fix pot_id variable
second@sam_data$pot_id <- sub("p","", second@sam_data$pot_id) 
second@sam_data$pot_id <- as.numeric(second@sam_data$pot_id) # this will let us sort by pot_id
second.sam <- data.frame(second@sam_data) %>%  # using arrange() in dplyr
  arrange(pot_id) %>%
  mutate(pot_id=paste0("p",pot_id)); head(second.sam); tail(second.sam)

# now make an identical dataframe to get variables from: 
sam <- data.frame(twice@sam_data) 
sam <- sam %>%
  filter(grepl(pattern = "_2", 
               sam$sample_id))

# and reassign variables: 
second.sam$host <- sam$host
second.sam$container <- sam$container
second.sam$interval <- sam$interval
second.sam$turn <- "second"; # add a variable to distinguish between merged and first samples

# convert til sample_data() object:
second.sam <- sample_data(second.sam); sample_names(second.sam)

# and reassign to twice.m phyloseq object: 
sample_data(second) <- second.sam; twice.m; head(second@sam_data)

# now we can merge the two phyloseq objects in order to normalize them:
fiphy <- merge_phyloseq(first, second); fiphy

# now we have to separate the new object into two dataframes before we combine
# them to calculate and plot fold increase: 
first <- subset_samples(fiphy, turn=="first")
first.df <- data.frame(cbind(estimate_richness(first, 
                                               measures = "Observed"),
                             first@sam_data)) %>% 
  rename("obs_first"="Observed",
         "pot_id_first"="pot_id",
         "turn_first"="turn") %>% 
  mutate(pot_id=pot_id_first); head(first.df)

second <- subset_samples(fiphy, turn=="second")
second.df <- data.frame(cbind(estimate_richness(second, 
                                                measures = "Observed"),
                              second@sam_data)) %>% 
  rename("obs_second"="Observed",
         "pot_id_second"="pot_id",
         "turn_second"="turn") %>% 
  mutate(pot_id=pot_id_second); head(second.df)

fidf <- full_join(first.df, second.df, by = "pot_id") %>% 
  select("obs_first", "obs_second", 
         "pot_id", "pot_id_first", "pot_id_second", 
         "turn_first", "turn_second")

fidf <- fidf %>% 
  mutate(foldincrease = obs_second/obs_first)
ggplot(fidf, aes(x=reorder(pot_id, -foldincrease), 
                 y=foldincrease,
                 fill="orange1", color="black")) + 
  geom_hline(yintercept = c(1, 1.25, 1.5, 1.75, 
                            2, 2.25, 2.5), 
             linetype = 2, color = "grey50") +
  geom_col(width=0.8, linewidth = 0.75) +
  scale_fill_manual(values = "orange1") +
  scale_color_manual(values = "black") +
  theme_classic() +
  theme(legend.position = "none", 
        axis.ticks.x = element_blank()) + 
  ylab(" ") +
  xlab("Sampled pots") +
  scale_x_discrete(labels = paste0(rep(" ", times=20))) +
  labs(title = "C") +
  scale_y_continuous(breaks = seq(0, 3, 0.25), limits = c(0, 3)) +
  coord_cartesian(ylim = c(1, 2.5)) -> fi.f; fi.f

# additionally, we will investigate whether there is a correlation between 
# the richness in the second sample and fold increase for each pot:
cormod <- cor.test(fidf$foldincrease, fidf$obs_second, 
                  method = "pearson"); cormod
ggplot(fidf, aes(x=obs_second, 
                 y=foldincrease,
                 fill="orange1")) + 
  geom_hline(yintercept = c(1, 1.25, 1.5, 1.75, 
                            2, 2.25, 2.5), 
             linetype = 2, color = "grey50") +
  theme_classic() +
  theme(legend.position = "none", 
  ) +
  scale_color_manual(values = "black") + 
  scale_fill_manual(values = "orange1") +
  geom_point(size=4,
             stroke = 1.25,
             pch=21) + 
  geom_smooth(method = "lm", 
              color="black", 
              fill="grey90", 
              linetype = 1) +
  xlab("Richness in both samples") +
  ylab(" ") + 
  labs(title = "D") +
  scale_y_continuous(breaks = seq(1, 3.5, 0.25), limits = c(0, 2.5)) +
  coord_cartesian(ylim = c(1, 2.5)) -> fi.point.f; fi.point.f

ggarrange(fi.p, fi.point.p, fi.f, fi.point.f, NULL, 
          widths = c(1, 1, 1, 1, 0.1), 
          nrow = 1)



