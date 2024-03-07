### 6_effect_of_sampling_strategy_on_contaminant_species_richness
### 2. sample volume

# we want to investigate the effect of sampling effort on the observed species 
# richness of contaminant plants and fungi identified by metabarcoding. 

# we want to compare contaminant species detection between two sample volumes. 
# unfortunately, all our samples were approx. 25 grams each. we will therefore
# simulate a double soil volume by merging samples that came from the same pots 
# (i.e., 10 pots per container) PRIOR to normalizing reads, and compare the
# observed species richness to that of pots only sampled once. 

### PLANTS

# we're going to merge samples from the same pots. to do this, we first create a 
# new phyloseq object with only the samples we're going to merge: 
twice <- subset_samples(plants, interval == "double"); twice
twice.m <- merge_samples(twice, group = "pot_id"); twice.m 

# reassign sample data to twice.m from twice:
twice.m@sam_data$pot_id <- rownames(twice.m@sam_data) # fix pot_id variable
twice.m@sam_data$pot_id <- sub("p","",twice.m@sam_data$pot_id) 
twice.m@sam_data$pot_id <- as.numeric(twice.m@sam_data$pot_id) # this will let us sort by pot_id - 
twice.m.sam <- data.frame(twice.m@sam_data) %>%  # - using arrange() in dplyr
  arrange(pot_id) %>%
  mutate(pot_id=paste0("p", pot_id))

# now make an identical dataframe to get variables from: 
sam <- data.frame(twice@sam_data) 
sam <- sam %>%
  filter(sample_id!="1_2", 
         sample_id!="2_2", 
         sample_id!="3_2",
         sample_id!="4_2",
         sample_id!="5_2")

# and reassign variables: 
twice.m.sam$host <- sam$host
twice.m.sam$container <- sam$container
twice.m.sam$interval <- sam$interval

# convert til sample_data() object:
twice.m.sam <- sample_data(twice.m.sam); sample_names(twice.m.sam)

# and reassign to twice.m phyloseq object: 
sample_data(twice.m) <- twice.m.sam; twice.m; head(twice.m@sam_data)


# now we want to get the merged samples back in the original dataset. we first
# create a new dataset with only samples from pots sampled once: 
once <- subset_samples(plants, interval == "single"); once

# and now we can merge the two phyloseq objects:
plants2 <- merge_phyloseq(twice.m, once); plants2

sort(sample_sums(plants2))
p2rar <- rarefy_even_depth(plants2, 
                          rngseed = 123, 
                          replace = FALSE, 
                          trimOTUs = TRUE, 
                          sample.size = 21915); plants2; p2rar

# now we're ready to have a look at observed species richness:
# make a dataframe for ease of testing and plotting:
obs.df <- data.frame(cbind(estimate_richness(p2rar, 
                                            measures = "Observed"), 
                          sample_data(p2rar))); str(obs.df)

ggdensity(obs.df$Observed)
ggqqplot(obs.df$Observed)

hist(obs.df$Observed); shapiro.test(obs.df$Observed) # not statistically significantly different from a normal distribution...
var.test(Observed ~ interval, obs.df) # the variance in the two groups is not significantly different from each other
summary(glm(Observed ~ interval, 
            data = obs.df, 
            family = gaussian))
summary(glm(Observed ~ host, 
            data = obs.df))
summary(glm(Observed ~ container, 
            data = obs.df))

summary(glm(Observed ~ interval * host, 
            data = obs.df))

mod <- glm(Observed ~ interval * host, data = obs.df); summary(mod)

emmeans::emmeans(mod, specs = list(pairwise ~ interval * host), adjust = "Tukey")

library(effects)

plot(allEffects(mod))

# make a nice boxplot to go with the results:
vol.plot.p <- obs.df %>% 
  ggplot(aes(x=interval, y=Observed, 
             fill=host)) + 
  geom_boxplot(alpha = 1, 
               linewidth = 0.75
               ) + 
  scale_fill_manual(values = c("chartreuse1", 
                               "chartreuse4")) +
  theme_classic() +
  ylab("Observed species richness") +
  xlab("Soil sample volume") + 
  scale_x_discrete(labels = c("Double", "Single")) +
  labs(title = "Soil sample volume", 
       subtitle = "Plants", 
       fill = "Volume") +
  annotate(geom = "line", 
           x = c(0.8, 1.2), 
           y = 15+0.2) +
  annotate(geom = "line", 
           x = c(0.8, 1.8), 
           y = 16+0.2) +
  annotate(geom = "line", 
           x = c(0.8, 2.2), 
           y = 17+0.2) +
  annotate(geom = "line", 
           x = c(1.2, 1.8), 
           y = 18+0.2) +
  annotate(geom = "line", 
           x = c(1.2, 2.2), 
           y = 19+0.2) +
  annotate(geom = "line", 
           x = c(1.8, 2.2), 
           y = 20+0.2) + 
  annotate(geom = "text", 
           label = c("**","*","**"), 
           x = c(1, 1.3, 1.4), 
           y = c(15.3,16.3,17.3), 
           size = 6) + 
  annotate(geom = "text", 
           label = c("n.s.","n.s.","n.s."), 
           x = c(1.5, 1.7, 2.0), 
           y = c(18.7,19.7,20.7), 
           size = 4) +
  theme(legend.title = element_blank(), 
        legend.position = c(0.85, 0.65),
        legend.background = element_rect(fill="white",
                                         linewidth=0.35,
                                         linetype="solid", 
                                         colour ="black"),
        legend.key.height = unit(0.5, 'cm'), 
        legend.key.width = unit(0.5, 'cm')); vol.plot.p

### FUNGI

# we're going to merge samples from the same pots. to do this, we first create a 
# new phyloseq object with only the samples we're going to merge: 
twice <- subset_samples(fungi, interval == "double"); twice
twice.m <- merge_samples(twice, group = "pot_id"); twice.m 

# reassign sample data in twice.m:
twice.m@sam_data$pot_id = rownames(twice.m@sam_data) # fix pot_id variable
twice.m@sam_data$pot_id = sub("p","",twice.m@sam_data$pot_id) 
twice.m@sam_data$pot_id = as.numeric(twice.m@sam_data$pot_id) # this will let us sort by pot_id
twice.m.sam = data.frame(twice.m@sam_data) %>%  # using arrange() in dplyr
  arrange(pot_id) %>%
  mutate(pot_id=paste0("p",pot_id))

# now make an identical dataframe to get variables from: 
sam <- data.frame(twice@sam_data) 
sam <- sam %>%
  filter(sample_id!="1_2", 
         sample_id!="2_2", 
         sample_id!="3_2",
         sample_id!="4_2",
         sample_id!="5_2")

# and reassign variables: 
twice.m.sam$host <- sam$host
twice.m.sam$container <- sam$container
twice.m.sam$interval <- sam$interval

# convert til sample_data() object:
twice.m.sam <- sample_data(twice.m.sam); sample_names(twice.m.sam)

# and reassign to twice.m phyloseq object: 
sample_data(twice.m) <- twice.m.sam; twice.m; head(twice.m@sam_data)


# now we want to get the merged samples back in the original dataset. we first
# create a new dataset with only samples from pots sampled once: 
once <- subset_samples(fungi, interval == "single"); once

# and now we can merge the two phyloseq objects:
fungi2 <- merge_phyloseq(twice.m, once); fungi2

sort(sample_sums(fungi2))
f2rar <- rarefy_even_depth(fungi2, 
                           rngseed = 123, 
                           replace = FALSE, 
                           trimOTUs = TRUE, 
                           sample.size = min(sample_sums(fungi2))); fungi2; f2rar

# now we're ready to have a look at observed species richness:
# make a dataframe for ease of testing and plotting:
obs.df <- data.frame(cbind(estimate_richness(f2rar, 
                                            measures = "Observed"), 
                          sample_data(f2rar))); str(obs.df)

hist(obs.df$Observed); shapiro.test(obs.df$Observed) # not statistically significantly different from a normal distribution...
var.test(log(Observed) ~ interval, obs.df) # the variance in the two groups is not significantly different from each other
summary(lm(Observed ~ interval, 
         data = obs.df))

mod <- lm(Observed ~ interval, data = obs.df); summary(mod)

emmeans::emmeans(mod, specs = list(pairwise ~ interval), adjust = "Tukey")

ggplot(obs.df, aes(x=interval, y=Observed, fill=interval)) + 
  geom_boxplot()


# make a nice boxplot to go with the results:
vol.plot.f <- obs.df %>% 
  ggplot(aes(x=interval, y=Observed, 
             fill=interval)) + 
  geom_boxplot(alpha = 1, 
               linewidth = 0.75) + 
  scale_fill_manual(values = c("orange1", 
                               "orange4")) +
  theme_classic() +
  ylab("Observed species richness") +
  xlab("Soil sample volume") + 
  scale_x_discrete(labels = c("Double", "Single")) +
  labs(title = " ", 
       subtitle = "Fungi", 
       fill = "Volume") +
  annotate(geom = "line", 
           x = c(1, 2), 
           y = 960) + 
  annotate(geom = "text", 
           label = c("***"), 
           x = c(1.5), 
           y = c(965), 
           size = 6) +
  theme(legend.title = element_blank(), 
        legend.position = "none"); vol.plot.f

ggarrange(vol.plot.p, 
          NULL,
          vol.plot.f, 
          widths = c(1, 0.2, 1),
          nrow = 1)
