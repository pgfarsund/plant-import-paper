### 5_effect_of_sampling_strategy_on_contaminant_species_richness
### 1. sample interval

# we want to investigate the effect of sampling effort on the observed species 
# richness of contaminant plants and fungi identified by metabarcoding. we'll do
# this in two steps:

# we want to compare contaminant species detection between two sampling 
# strategies for obtaining 20 samples per container: sampling 20 pots once, or
# sampling 10 pots twice. 
# we're interested in 1) the total number of species identified by both 
# strategies, and 2) estimates of total richness based on the two strategies:

library(iNEXT)
sort(sample_sums(plants))

# normalise:
p2rar <- rarefy_even_depth(plants, 
                          rngseed = 123, 
                          replace = FALSE, 
                          trimOTUs = TRUE, 
                          sample.size = 21915); plants; p2rar

# identify pots that were samples twice but lost a sample during normalizing: 
ggplot(data.frame(p2rar@sam_data), 
       aes(x=pot_id)) + 
  geom_bar() + 
  facet_wrap(~interval, scales = "free_x")

# pots 14, 20, 46, and 48 lost one sample. filter them out of the dataset:
p2rar <- subset_samples(p2rar, pot_id!="p17"); p2rar
p2rar <- subset_samples(p2rar, pot_id!="p20"); p2rar
p2rar <- subset_samples(p2rar, pot_id!="p46"); p2rar
p2rar <- subset_samples(p2rar, pot_id!="p48"); p2rar

table(p2rar@sam_data$interval)

# to run iNEXT, we need a list of objects to interpolate/extrapolate richness in. 
# we want one curve for each sampling effort per container, i.e. four curves in total. 
# start by sub-setting containers and sampling efforts: 

cont1 <- subset_samples(p2rar, container == "1") # subset
cont1 <- prune_taxa(taxa_sums(cont1)>0, cont1); cont1 # prune taxa with 0 abundance

# make an asv-by-sample data frame with presence/absence (1/0) data: 
cont1inext <- data.frame(t(cont1@otu_table)) %>% 
  mutate_if(is.numeric, ~1 *(. !=0)) # transform values >1 into 1

# make similar asv-by-sample data frames per container and sampling interval: 
c1s <- subset_samples(cont1, interval=="single"); c1s <- prune_taxa(taxa_sums(c1s)>0, c1s); c1s
c1s <- data.frame(t(c1s@otu_table)) %>% 
  mutate_if(is.numeric, ~1 *(. !=0))
c1d <- subset_samples(cont1, interval=="double"); c1d <- prune_taxa(taxa_sums(c1d)>0, c1d); c1d
c1d <- data.frame(t(c1d@otu_table)) %>% 
  mutate_if(is.numeric, ~1 *(. !=0))

cont2 = subset_samples(p2rar, container == "2"); cont2 <- prune_taxa(taxa_sums(cont2)>0, cont2); cont2
cont2inext = data.frame(t(cont2@otu_table)) %>% 
  mutate_if(is.numeric, ~1 *(. !=0))

c2s <- subset_samples(cont2, interval=="single"); c2s <- prune_taxa(taxa_sums(c2s)>0, c2s); c2s
c2s <- data.frame(t(c2s@otu_table)) %>% 
  mutate_if(is.numeric, ~1 *(. !=0))
c2d <- subset_samples(cont2, interval=="double"); c2d <- prune_taxa(taxa_sums(c2d)>0, c2d); c2d
c2d <- data.frame(t(c2d@otu_table)) %>% 
  mutate_if(is.numeric, ~1 *(. !=0))

# make a list compatible with iNEXT:
conts <- list(c1s, c1d, 
              c2s, c2d)

# run iNEXT:
conts.out.p <- iNEXT::iNEXT(conts, 
                            q = 0, 
                            datatype = "incidence_raw", 
                            knots = 40, 
                            endpoint = 40,
                            se = TRUE)
inext.res.p <- conts.out.p

# fix the names of our assemblages to match our dataset: 
inext.res.p$DataInfo$Assemblage <- sub("assemblage1", "Cont.1 one/pot", inext.res.p$DataInfo$Assemblage)
inext.res.p$iNextEst$size_based$Assemblage <- sub("assemblage1", "Cont.1 one/pot", inext.res.p$iNextEst$size_based$Assemblage)
inext.res.p$iNextEst$coverage_based$Assemblage <- sub("assemblage1", "Cont.1 one/pot", inext.res.p$iNextEst$coverage_based$Assemblage)
inext.res.p$AsyEst$Assemblage <- sub("assemblage1", "Cont.1 one/pot", inext.res.p$AsyEst$Assemblage)

inext.res.p$DataInfo$Assemblage <- sub("assemblage2", "Cont.1 two/pot", inext.res.p$DataInfo$Assemblage)
inext.res.p$iNextEst$size_based$Assemblage <- sub("assemblage2", "Cont.1 two/pot", inext.res.p$iNextEst$size_based$Assemblage)
inext.res.p$iNextEst$coverage_based$Assemblage <- sub("assemblage2", "Cont.1 two/pot", inext.res.p$iNextEst$coverage_based$Assemblage)
inext.res.p$AsyEst$Assemblage <- sub("assemblage2", "Cont.1 two/pot", inext.res.p$AsyEst$Assemblage)

inext.res.p$DataInfo$Assemblage <- sub("assemblage3", "Cont.2 one/pot", inext.res.p$DataInfo$Assemblage)
inext.res.p$iNextEst$size_based$Assemblage <- sub("assemblage3", "Cont.2 one/pot", inext.res.p$iNextEst$size_based$Assemblage)
inext.res.p$iNextEst$coverage_based$Assemblage <- sub("assemblage3", "Cont.2 one/pot", inext.res.p$iNextEst$coverage_based$Assemblage)
inext.res.p$AsyEst$Assemblage <- sub("assemblage3", "Cont.2 one/pot", inext.res.p$AsyEst$Assemblage)

inext.res.p$DataInfo$Assemblage <- sub("assemblage4", "Cont.2 two/pot", inext.res.p$DataInfo$Assemblage)
inext.res.p$iNextEst$size_based$Assemblage <- sub("assemblage4", "Cont.2 two/pot", inext.res.p$iNextEst$size_based$Assemblage)
inext.res.p$iNextEst$coverage_based$Assemblage <- sub("assemblage4", "Cont.2 two/pot", inext.res.p$iNextEst$coverage_based$Assemblage)
inext.res.p$AsyEst$Assemblage <- sub("assemblage4", "Cont.2 two/pot", inext.res.p$AsyEst$Assemblage)
inext.res.p$DataInfo$container <- c("1","1","2","2")

# finally we plot our curves: 
p3 <- ggiNEXT(inext.res.p, type = 1, se = F) +
  #geom_borderline(size=1, bordercolor = "black", alpha = 1, linetype = 1) +
  theme_classic() +
  scale_color_manual(values = c("chartreuse1", "chartreuse2", "chartreuse3", "chartreuse4"), 
                     labels = c("\nContainer 1\n20 pots\n", 
                                "\nContainer 1\n10 pots\n", 
                                "\nContainer 2\n20 pots\n", 
                                "\nContainer 2\n10 pots\n")) +
  scale_fill_manual(values = c("chartreuse1", "chartreuse2", "chartreuse3", "chartreuse4"), 
                     labels = c("\nContainer 1\n20 pots\n", 
                                "\nContainer 1\n10 pots\n", 
                                "\nContainer 2\n20 pots\n", 
                                "\nContainer 2\n10 pots\n")) +
  scale_shape_manual(values = c(21, 22, 23, 24),
                     labels = c("\nContainer 1\n20 pots\n", 
                                "\nContainer 1\n10 pots\n", 
                                "\nContainer 2\n20 pots\n", 
                                "\nContainer 2\n10 pots\n")) +
  theme(#legend.position = c(0.7, 0.2), 
    legend.position = "right", 
    legend.direction = "vertical", 
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.title.x = element_text(),
    axis.title.y = element_text(), 
    legend.text = element_text(size = 10)) +
  ylab("Species richness") +
  labs(title = " ", 
       subtitle = "A.") + 
  xlab("Number of samples") +
  guides(color = guide_legend(direction = "vertical",
                              title = NULL, 
                              nrow = 4, byrow = T),
         fill = guide_legend(direction = "vertical",
                             title = NULL, 
                             nrow = 4, byrow = T),
         shape = guide_legend(direction = "vertical",
                              title = NULL, 
                              nrow = 4, byrow = T),
         line = guide_legend(direction = "vertical", 
                             title = NULL, 
                             nrow = 4, byrow = T)); p3
psac <- ggplot_build(p3)

psac$data[[1]]$size <- 5
psac$data[[1]]$stroke <- c(1,1,1,1)
psac$data[[1]]$colour <- c("black","black","black","black")
psac$data[[1]]$alpha <- c(1, 1, 1, 1)
psac$data[[1]]$fill <- c("chartreuse1", "chartreuse2", "chartreuse3", "chartreuse4")
psacplot <- ggplot_gtable(psac); plot(psacplot)

library(ggplotify)
psacplot <- as.ggplot(psac); psacplot # pdf landscape 4 x 5

#ggarrange(psac, volplants, nrow = 1)
#
results <- data.frame(inext.res.p$AsyEst) %>% 
  filter(Diversity == "Species richness") %>% 
  dplyr::select(-Diversity) %>% 
  mutate(across(c("Estimator", 
                  "s.e.",
                  "LCL",
                  "UCL"), round, 0)) %>% 
  print()

p3 + annotate(geom = "text", x = 10, y = 15, hjust = 0, size = 3,
              label = insight::export_table(results)) -> p3; p3
