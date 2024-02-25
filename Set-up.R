################################################################################
#                          Required Libraries                                  #
{
library(tidyr)
library(readr)
library(ggplot2)
library(gganimate)
library(dplyr)
library(ggpubr)
library(plotly)
library(tidyverse)
library(ggtext)
library(RColorBrewer)
library(vegan)
library(readxl)
}
################################################################################
#Reading in the final.opti_mcc.shared file from mothur output           
################################################################################
#first let's read in the metadata file you have on the samples
metadata <- metadata
#read in the final.opti_mcc.shared tsv file from mothur output and making it a tibble that has count of each Otu for each sample
otu_counts <- read_tsv("/Users/trentonwolfe/UofH-Collab/IBZ Manuscript/mothur/mothur.out/final.opti_mcc.shared")%>%
  select(-label, -numOtus)%>%
  pivot_longer(-Group, names_to = "otu", values_to = "count")
otu_counts #print the tibble

#extract the group names from the otu_counts vector so I can ensure my metadata only has these and work in the external data for these mice
groups <- as.data.frame(unique(otu_counts$Group))
write_excel_csv(groups, file="groups.csv")

#now, let's read in the final.opti_mcc.0.03.cons.taxonomy file from mothur to join counts and taxonomy together for rel abund
taxonomy<- read_tsv("/Users/trentonwolfe/UofH-Collab/IBZ Manuscript/mothur/mothur.out/final.opti_mcc.0.03.cons.taxonomy")%>%
  rename_all(tolower)%>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)",""),
         axonomy = str_replace(taxonomy, ";$", ""))%>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"), 
           sep=";")

metadata <- metadata
#cool, now combine otu_counts with your metadata for each  the taxonomy for the OTUs that we just did above and calculate a relative abundance for each OTU
otu_rel_abund <- inner_join(otu_counts, metadata, by="Group")%>%
  inner_join(., taxonomy, by="otu")%>%
  group_by(Group)%>%
  mutate(rel_abund = count / sum(count))%>%
  ungroup()%>%
  select(-count)%>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus", "otu"),
               names_to="level",
               values_to="taxon")%>%
  mutate(taxon = gsub("\\s*\\(\\d+\\)$", "", taxon))
#this otu_rel_abund df will be used in the Rel-Abund file

#now let's read in the metrics I calculated in mothur 
metrics <- read_tsv("/Users/trentonwolfe/Downloads/final.opti_mcc.groups.summary")%>%
  rename(Group = group)

#now I'm going to combine in the diversity metrics I calculated in mothur and paste them in with the metadata file
metadata2 <- inner_join(metrics,metadata)
#this "metadata2" file can now be used to do alpha diversity analysis as it now has all the mouse metadata and the alpha diversity metrics that were calculated in mothur for the mice



#make a color vector to use for all future plots
unique(metadata2$`abx-grp`)
abx.colors <- c("Baseline 1" = "grey", 
              "Baseline 2" = "darkgrey",
              "NDC" = "black", 
              "IBZ" = "#F267B1",
              "FDX" = "#1C8F63",
              "VAN" = "#0047AB",
              "MTZ" = "#946317")

abx.col <- c("NDC" = "black",
             "FDX" = "#1C8F63",
             "IBZ" = "#F267B1",
             "VAN" = "#0047AB",
             "MTZ" = "#946317")





