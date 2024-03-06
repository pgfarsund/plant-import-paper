# 4_contaminant_species_diversity

# we want to investigate 

# normalise data: 
p.rar <- rarefy_even_depth(plants, 
                           rngseed = 123, 
                           replace = FALSE, 
                           trimOTUs = TRUE,
                           sample.size = 21915); plants; p.rar

# create a data.frame with observed richness and metadata:
left_join(data.frame(sample_data(p.rar), 
                     sam_id = sample_names(p.rar)),
          data.frame(estimate_richness(p.rar, measures = "Observed"), 
                     sam_id = sample_names(p.rar))) %>% 
  column_to_rownames(var = "sam_id") -> obs

# we do some plotting to explore the data:
obs %>% 
  ggplot(aes(x=container, y=Observed, 
             fill=host)) + 
  geom_boxplot()

# it seems like there is a difference in observed species richness between thuja
# and taxus in container 1 but not in container 2. this indicates an interaction 
# effect between container and host. 

# before fitting our data to a model, we need to figure out what distribution
# best fits our data. we use the package fitdistrplus for this:
library(fitdistrplus) # load package

# we found it difficult to fit the data to a suitable distribution. by square
# root transforming the data, they're not significantly different from a normal
# distribution:
plotdist(sqrt(obs$Observed), 
         histo = TRUE, demp = TRUE) # intial plotting

descdist(sqrt(obs$Observed), 
         boot = 1000) # skewness-kurtosis plot

x=fitdist(sqrt(obs$Observed), 
          distr = "norm", 
          method = "mle"); summary(x); plot(x) 

gofstat(x) # not statistically significant different from a normal distribution. 

# we fit a linear model to test the effect of container and host plant on the
# on observed species richness. include pot_id as a random effect to control
# for the fact that ca. half of the pots were sampled twice: 
library(lmerTest) # for fitting the model
library(effects) # for plotting effects of predictor variables

p.mod = lmerTest::lmer(formula = sqrt(Observed) ~ container * host + (1|pot_id),
                       data = obs); summary(p.mod); plot(allEffects(p.mod))

# we'll use DHARMa to validate the model's residuals:
library(DHARMa)
testDispersion(p.mod)

sim.res.p.mod <- simulateResiduals(p.mod, 
                                  n = 1000,
                                  plot = TRUE) # these simulated residuals look good.

# review the model results:
p.mod.res <- summary(p.mod); p.mod.res

# host and container:host had a significant effect on plant richness in our samples. 
# perform a post-hoc test with emmeans package:
library(multcomp)
ph <- glht(p.mod, linfct = mcp(list(container * host = "Tukey"))); summary(ph)

library(emmeans)
p.posthoc <- emmeans(p.mod, 
                     list(pairwise ~ container*host), 
                     adjust = "Tukey"); p.posthoc # two significant pairwise comparisons

data.frame(p.posthoc$`pairwise differences of container, host`) %>%
  filter(p.value < 0.05) 

pmodbox <- obs %>% 
  ggplot(aes(x=container, y=Observed, 
             fill=host)) + 
  geom_boxplot() + 
  geom_boxplot(alpha = 0.6) + 
  scale_fill_manual(values = c("darkolivegreen1", 
                               "darkolivegreen")) +
  theme_classic() +
  ylab("Observed species richness") +
  scale_x_discrete(labels = c("Container 1", "Container 2")) +
  labs(title = "A.") +
  annotate(geom = "line", 
           x = c(0.8, 1.2), 
           y = 20+0.2) +
  annotate(geom = "line", 
           x = c(0.8, 1.8), 
           y = 19+0.2) +
  annotate(geom = "line", 
           x = c(0.8, 2.2), 
           y = 18+0.2) +
  annotate(geom = "line", 
           x = c(1.2, 1.8), 
           y = 17+0.2) +
  annotate(geom = "line", 
           x = c(1.2, 2.2), 
           y = 16+0.2) +
  annotate(geom = "line", 
           x = c(1.8, 2.2), 
           y = 15+0.2) + 
  annotate(geom = "label", 
           label = c("***","n.s.","n.s."), 
           x = c(1, 1.3, 1.4), 
           y = c(20,19,18), 
           size = 4) + 
  annotate(geom = "label", 
           label = c("**","*","n.s."), 
           x = c(1.5, 1.7, 2.0), 
           y = c(17,16,15), 
           size = 4, margin) + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.85, 0.55),
        legend.background = element_rect(fill="grey95",
                                         linewidth=0.35,
                                         linetype="solid", 
                                         colour ="black"),
        axis.title.x = element_blank(),
        legend.key.height = unit(0.5, 'cm'), 
        legend.key.width = unit(0.5, 'cm')); pmodbox
#