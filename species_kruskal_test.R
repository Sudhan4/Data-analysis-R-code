rm(list = ls())

library(readxl); library(changepoint); library(rstatix); library(reshape); library(tidyverse); 
library(ggpubr); library(datarium); library(dbplyr); library(dplyr); library(tidyr);
library("ggsci"); library(RColorBrewer); library(ggplot2); library(ggsignif); library(devtools);
library(pairwiseAdonis); library(ggplot2); library(grid); library(ggforce); library(microbiome);
library(dplyr); library(vegan); library(lsr); library(effsize)
library(broom)

taxonomy_rel_abun <- as.data.frame(read.csv("taxonomy_features.csv", row.names = 1)) 

rownames(taxonomy_rel_abun) <- gsub("^.*__","",rownames(taxonomy_rel_abun))
rownames(taxonomy_rel_abun) <- gsub("_"," ",rownames(taxonomy_rel_abun))


metadata <- as.data.frame(read.table("metadata_new.txt")) 

metadata$group <- as.factor(metadata$group) 

transp_taxonomy_rel_abun <- as.data.frame(t(taxonomy_rel_abun)) %>% 
  select(-starts_with("GGB")) #Exclude uncharacterized species

species_metadata <- as.data.frame(merge(transp_taxonomy_rel_abun,metadata, by = 0)) 


composite <- species_metadata %>% 
  pivot_longer(2:126, 
               names_to = "species", 
               values_to = "rel_abun") %>% 
  select(species,rel_abun, time, group, strain)

significant_species <- composite %>% nest(data = -species) %>% 
  mutate(test = map(.x = data, ~kruskal.test(rel_abun~group, data=.x) %>% tidy)) %>% 
  unnest(test) %>% 
  filter(p.value<0.05)
