# we will make two species-by-pot matrices with presence/absence (1/0) data, one 
# for the first samples in every pot and one for the second pot. we'll then merge
# the two matrices, leaving every cell with 0, 1, or 2 identifications per species
# per pot. for every pot, each species will then be
#
#       0 = not identified 
#       1 = identified in one sample
#       2 = identified in both samples
# 
# this will make it easy for us to examine the recapture rate of the different species. 

twice = subset_samples(plants, interval == "double"); twice # subset pots sampled twice

first = subset_samples(twice, grepl(pattern = "_1", # subset first samples
                                    x = twice@sam_data$sample_id)); first = prune_taxa(taxa_sums(first)>0,
                                                                                       first); first
first@sam_data$turn = "first"; head(first@sam_data) # assign variable to distinguish from "second" samples
first@sam_data$row = paste0(first@sam_data$pot_id, # pot ID
                            "_",                   # _
                            first@sam_data$turn)   # turn

second = subset_samples(twice, grepl(pattern = "_2", # second samples
                                     x = twice@sam_data$sample_id)); second = prune_taxa(taxa_sums(second)>0,
                                                                                         second); second
second@sam_data$turn = "second"; head(second@sam_data)
second@sam_data$row = paste0(second@sam_data$pot_id, # pot ID
                             "_",                    # _
                             second@sam_data$turn)   # turn

twice = merge_phyloseq(first, 
                       second); twice # merge back together

sample_names(twice) = twice@sam_data$row; head(twice@sam_data) # assign the new sample names

colnames(twice@otu_table) = twice@tax_table[,6]

# now we make the new data frame with presence/absence data:
pres.abs.df = data.frame(twice@otu_table) %>% 
  mutate_if(is.numeric, ~1 * (. > 0)) 

pres.abs.df.t = data.frame(t(pres.abs.df)) # transpose to make pots columns

# split first and second samples and clean up the colnames before summarizing:
first = data.frame(pres.abs.df.t) %>% 
  dplyr::select(matches("_first"))
colnames(first) = sub("_first", "", colnames(first))

second = data.frame(pres.abs.df.t) %>% 
  dplyr::select(matches("_second"))
colnames(second) = sub("_second", "", colnames(second))

# by summarizing the two data frames, we'll be able to see if species were identified
# 0, 1, or 2 times per pot:
pres.abs.df.fin = data.frame(first+second)

# this is where you find the overall recapture rate of 70.3 % --- that does not add up....
# the recapture rate has to be the number of 2 divided by the number of 1s + 2s
table(pres.abs.df.fin==0) # 175 identifications in total (i.e. != 0)
table(pres.abs.df.fin==1) # species X was identified in only one sample per pot a total of 123 times
table(pres.abs.df.fin==2) # species X was identified in both samples per pot a total of 52 times

# check recapture rate for two known alien species to Norway:
con = pres.abs.df.fin["Erigeron_canadensis",] # E. canadensis is synonymous to Conyza canadensis
table(con==0) 
table(con==1) 
table(con==2)

acer = pres.abs.df.fin[rownames(pres.abs.df.fin) %in% "Acer_pseudoplatanus", ]
acer = pres.abs.df.fin["Acer_pseudoplatanus",]
table(acer==0) 
table(acer==1) 
table(acer==2)










# we will make two species-by-pot matrices with presence/absence (1/0) data, one 
# for the first samples in every pot and one for the second pot. we'll then merge
# the two matrices, leaving every cell with 0, 1, or 2 identifications per species
# per pot. for every pot, each species will then be
#
#       0 = not identified 
#       1 = identified in one sample
#       2 = identified in both samples
# 
# this will make it easy for us to examine the recapture rate of the different species. 

twice = subset_samples(fungi, interval == "double"); twice # subset pots sampled twice

first = subset_samples(twice, grepl(pattern = "_1", # subset first samples
                                    x = twice@sam_data$sample_id)); first = prune_taxa(taxa_sums(first)>0,
                                                                                       first); first
first@sam_data$turn = "first"; head(first@sam_data) # assign variable to distinguish from "second" samples
first@sam_data$row = paste0(first@sam_data$pot_id, # pot ID
                            "_",                   # _
                            first@sam_data$turn)   # turn

second = subset_samples(twice, grepl(pattern = "_2", # second samples
                                     x = twice@sam_data$sample_id)); second = prune_taxa(taxa_sums(second)>0,
                                                                                         second); second
second@sam_data$turn = "second"; head(second@sam_data)
second@sam_data$row = paste0(second@sam_data$pot_id, # pot ID
                             "_",                    # _
                             second@sam_data$turn)   # turn

twice = merge_phyloseq(first, 
                       second); twice # merge back together

sample_names(twice) = twice@sam_data$row; head(twice@sam_data) # assign the new sample names

colnames(twice@otu_table) = twice@tax_table[,7]

# now we make the new data frame with presence/absence data:
pres.abs.df = data.frame(twice@otu_table) %>% 
  mutate_if(is.numeric, ~1 * (. > 0)) 

pres.abs.df.t = data.frame(t(pres.abs.df)) # transpose to make pots columns

# split first and second samples and clean up the colnames before summarizing:
first = data.frame(pres.abs.df.t) %>% 
  dplyr::select(matches("_first")) # friggin' MASS...
colnames(first) = sub("_first", "", colnames(first))

second = data.frame(pres.abs.df.t) %>% 
  dplyr::select(matches("_second"))
colnames(second) = sub("_second", "", colnames(second))

# by summarizing the two data frames, we'll be able to see if species were identified
# 0, 1, or 2 times per pot:
pres.abs.df.fin = data.frame(first+second)

# this is where you find the overall recapture rate of 70.3 % 
table(pres.abs.df.fin==0) # 175 identifications in total (i.e. != 0)
table(pres.abs.df.fin==1) # species X was identified in only one sample per pot a total of 123 times
table(pres.abs.df.fin==2) # species X was identified in both samples per pot a total of 52 times

# check recapture rate for two known alien species to Norway:

crypt = pres.abs.df.fin["Cryptostroma_corticale",] # subset C. corticale
table(crypt==0) # C. corticale was unidentified in 13 pots
table(crypt==1) # C. corticale was identified in one sample in 6 pots
table(crypt==2) # C. corticale was identified twice in one pot

mut = pres.abs.df.fin["Mutinus_ravenelii",] # subset M. ravenellii
table(mut==0) 
table(mut==1) 
table(mut==2)


