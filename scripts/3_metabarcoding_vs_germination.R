### 3_metabarcoding_vs_germination

# we identified 33 plant species with metabarcoding and 30 with germination. five 
# species were identified by both methods, adding up to 58 contaminant plant species
# identified all together. we also checked whether identified plants were registered 
# as products in the online assortments of Noviflora and Plantasjen in the fall of 
# 2021 (the year of sampling). we will make two stacked bar charts to show how many 
# species fall into each category, while also showing how many species were identified
# with metabarcoding, germination, or both: 

df <- read.csv(here("data", "fig_2_dataset.csv"),
               sep = ";"); str(df) # change path 

my_colors = c("#CC79A7",
              "#009E73",
              "#E69F00")

a <- df %>% 
  mutate(prev_det=sub("TRUE", "Yes", prev_det)) %>% 
  mutate(prev_det=sub("FALSE", "No", prev_det)) %>% 
  ggplot(aes(x=factor(prev_det, 
                      levels = c("Yes", 
                                 "No")),
             color=factor(id_method,
                          levels = c("Metabarcoding",
                                     "Germination",
                                     "Both")),
             fill=factor(id_method,
                         levels = c("Metabarcoding",
                                    "Germination",
                                    "Both")))) +
  geom_hline(yintercept = c(10, 20, 30, 40), 
             linetype = 2, 
             linewidth = 0.5) +
  geom_bar(position = "stack",
           width = 0.55, 
           linewidth = 0.75, 
           alpha = 1) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = rep("black", times=3)) +
  theme_classic(#base_size = 15
                ) + 
  xlab("Previously identified?") + 
  ylab("Number of plant species") + 
  labs(fill = "Identification\nmethod", 
       color = "Identification\nmethod", 
       subtitle = "a.") +
  theme(legend.title = element_text(#size = 13
                                    )) +
  scale_y_continuous(breaks = seq(0, 45, 10), 
                     limits = c(0, 45)); a

b <- df %>% 
  mutate(ornamental=sub("TRUE", "Yes", ornamental)) %>% 
  mutate(ornamental=sub("FALSE", "No", ornamental)) %>% 
  ggplot(aes(x=factor(ornamental, 
                      levels = c("Yes", 
                                 "No")),
             color=factor(id_method,
                          levels = c("Metabarcoding",
                                     "Germination",
                                     "Both")),
             fill=factor(id_method,
                         levels = c("Metabarcoding",
                                    "Germination",
                                    "Both")))) +
  geom_hline(yintercept = c(10, 20, 30, 40), 
             linetype = 2, 
             linewidth = 0.5) +
  geom_bar(position = "stack",
           width = 0.55, 
           linewidth = 0.75, 
           alpha = 1) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = rep("black", times=3)) +
  theme_classic(#base_size = 15
                ) + 
  xlab("Ornamental?") + 
  ylab(" ") + 
  labs(fill = "Identification\nmethod", 
       color = "Identification\nmethod", 
       subtitle = "b.") +
  theme(legend.title = element_text(#size = 13
                                    )) +
  scale_y_continuous(breaks = seq(0, 55, 10), 
                     limits = c(0, 45)); b

fig.2 <- ggarrange(a, b, 
                   common.legend = TRUE,
                   legend = "right", 
                   nrow = 1); fig.2 # 3.0 x 5.0
