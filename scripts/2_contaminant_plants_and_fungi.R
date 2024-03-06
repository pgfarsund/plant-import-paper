### 2_contaminant_plants_and_fungi

library(phyloseq)
library(tidyverse)
library(ggpubr)

# we identified many contaminant plants and fungi that were native, known aliens, 
# and potential alien species to Norway. we've grouped these species as "Native", 
# "Alien", "Swe., Fin., Den." (i.e. native to Norway's neighboring countries), 
# (native to) "Europe", and (native to) "Non-Europe". we will make a stacked bar
# chart to show the relative frequencies of these different categories:

# construct dataframe: 
df = data.frame("freq" = c(30.30, 45.45, 0, 12.12, 12.12, 0, # start with plant frequencies
                           39.53, 0.81, 22.32, 22.89, 10.88, 3.57), # then fungal frequencies
                "group" = rep(c("Plants",
                                "Fungi"), each=6),
                "cats" = c(rep(c("Native",
                                 "Alien",
                                 "Swe., Fin., Den.",
                                 "Europe", 
                                 "Non-Europe",
                                 "NA"), times=2)))

# make Figure 1 in the manuscript:
fig.1 = df %>%
  ggplot(aes(x=factor(group, 
                      level=c("Plants",
                              "Fungi")), # order the groups on the x-axis
             y=freq, 
             fill=factor(cats, 
                         levels = c("Native", # we want the stacked bars to have this specific order
                                    "Alien", 
                                    "Swe., Fin., Den.",
                                    "Europe", 
                                    "Non-Europe", 
                                    "NA")),
             color=factor(cats, 
                          levels = c("Native", 
                                     "Alien", 
                                     "Swe., Fin., Den.",
                                     "Europe", 
                                     "Non-Europe", 
                                     "NA")))) + 
  geom_col(position = "fill", 
           width = 0.55) + 
  theme_classic(base_size = 15) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", # some nice colors to fill the bars
                               "#F0E442", "#CC79A7", "#999999")) +
  scale_color_manual(values = rep("black", times = 6)) +
  theme(legend.position = "right") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank()) +
  annotate(geom = "text", 
           label = c("30.3 %", 
                     "45.5 %",
                     "12.1 %",
                     "12.1 %", 
                     "39.5 %",
                     "0.8 %",
                     "22.3 %", 
                     "22.9 %", 
                     "10.9 %",
                     "3.6 %"),
           x = rep(c(1,
                     2,
                     2.45,
                     2,
                     2.45), # we want "0.8 %" and "3.6 %" outside the actual bar, which makes this a bit tricky
                   times = c(4,
                             1,
                             1,
                             3,
                             1)), 
           y = c(0.84, 0.45, 0.18, 0.06, # plant percent labels
                 0.8, 0.6, 0.5, 0.26, 0.1, 0.02), # fungi percent labels
           size = 4, color = "black", fontface = 1) + 
  annotate(geom = "segment", 
           x = rep(2.25,
                   times=2), 
           xend = rep(2.315, 
                      times=2), 
           y = c(0.60,
                 0.02), 
           yend = c(0.60,
                    0.02)) + 
  scale_y_continuous(labels = scales::percent); fig.1 # save as .pdf landscape 5 x 6