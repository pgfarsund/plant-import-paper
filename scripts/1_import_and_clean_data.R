### 1_import_and_clean_data

# plants
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(here)

### plants
plants <- readRDS(here("data", "plants-phyloseq-object.RDS"))

# manipulate the sample data:
sam <- data.frame(plants@sam_data) %>% 
  dplyr::select("sample","Location_name",            # there is a lot of data we don't need.
                "Specific_host",
                "Sample_ID") %>% 
  rename(c("container"="Location_name",              # rename columns for simplicity.
           "host"="Specific_host", 
           "sample_id"="Sample_ID")) %>% 
  mutate(sample=sub("PLANTEIMP_21_", "", sample),    # easier to interpret sample number.
         sample=sub("_ITS_S2_ITS4", "", sample)) %>% 
  mutate(container = sub("Bil ", "", container)) %>% # remove "Bil " from the container ID.
  mutate(pot_id = c(rep(paste0("p", 1:5), each=2),   # add an ID to each pot - this will be
                    rep(paste0("p", 6:15), each=1),  # identical between samples taken from
                    rep(paste0("p", 16:20), each=2), # the same pot.
                    rep(paste0("p", 21:30), each=1),
                    rep(paste0("p", 31:35), each=2),
                    rep(paste0("p", 36:45), each=1),
                    rep(paste0("p", 46:50), each=2),
                    rep(paste0("p", 51:60), each=1))) %>%
  mutate(interval = c(rep("double", times=10),
                      rep("single", times=10),
                      rep("double", times=10),
                      rep("single", times=10),
                      rep("double", times=10),
                      rep("single", times=10), 
                      rep("double", times=10),
                      rep("single", times=10))); str(sam) 
sample_data(plants) <- sample_data(sam); rm(sam); plants

# for this study we're only interested in vascular plants. our data set contains
# some bryophytes and green algae, additionally we're not interested in Taxus sp. 
# or Thuja sp., because these were ornamentals in the sampled soil. 
# remove non-focal taxa:
plants <- subset_taxa(plants, phylum=="Streptophyta")
plants <- subset_taxa(plants, class!="Chlorophyceae")
plants <- subset_taxa(plants, order!="Splachnales")
plants <- subset_taxa(plants, order!="Bryales")
plants <- subset_taxa(plants, order!="Dicranales")
plants <- subset_taxa(plants, order!="Funariales")
plants <- subset_taxa(plants, order!="Hypnales")
plants <- subset_taxa(plants, order!="Klebsormidiaceae")
plants <- subset_taxa(plants, order!="Pottiales")
plants <- subset_taxa(plants, order!="Pseudoditrichales")
plants <- subset_taxa(plants, species!="Thuja_sp.")
plants <- subset_taxa(plants, species!="Thuja_occidentalis")
plants <- subset_taxa(plants, species!="Taxus_sp.")
plants <- subset_taxa(plants, species!="Taxus_canadensis"); plants



### fungi

# load data:
fungi <- readRDS(here("data", "fungi-phyloseq-object.RDS")) 
fungi <- subset_samples(fungi, negative_control=="FALSE") # remove negative controls

# manipulate the sample data:
sam = data.frame(fungi@sam_data) %>% 
  dplyr::select("sample","Location_name",            # there is a lot of data we don't need
                "Specific_host",
                "Sample_ID") %>%
  rename(c("container"="Location_name",              # rename columns for simplicity
           "host"="Specific_host", 
           "sample_id"="Sample_ID")) %>% 
  mutate(sample=sub("PLANTEIMP_21_", "", sample),    # easier to interpret sample number
         sample=sub("_fITS7_ITS4", "", sample)) %>% 
  mutate(container = sub("Bil ", "", container)) %>% # remove "Bil " from the container ID
  mutate(pot_id = c(rep(paste0("p", 1:5), each=2),   # add an ID to each pot - this will be
                    rep(paste0("p", 6:15), each=1),  # identical between samples taken from
                    rep(paste0("p", 16:20), each=2), # the same pot
                    rep(paste0("p", 21:30), each=1),
                    rep(paste0("p", 31:35), each=2),
                    rep(paste0("p", 36:45), each=1),
                    rep(paste0("p", 46:50), each=2),
                    rep(paste0("p", 51:60), each=1))) %>%
  mutate(interval = c(rep("double", times=10),
                      rep("single", times=10),
                      rep("double", times=10),
                      rep("single", times=10),
                      rep("double", times=10),
                      rep("single", times=10), 
                      rep("double", times=10),
                      rep("single", times=10)))
sample_data(fungi) <- sample_data(sam); rm(sam); fungi
