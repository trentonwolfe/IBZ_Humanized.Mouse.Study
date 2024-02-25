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
#Relative Abundance Analysis    
################################################################################
#To start, we need to generate a table of otu counts w/ metadata and tax. labs
otu_rel_abund <- inner_join(metadata, otu_counts, by="Group")%>%
  inner_join(., taxonomy, by="otu")%>%
  group_by(Group)%>%
  mutate(rel_abund = count / sum(count))%>%
  ungroup()%>%
  select(-count)%>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus", "otu"),
               names_to="level",
               values_to="taxon")%>%
  mutate(taxon = gsub("\\s*\\(\\d+\\)$", "", taxon))

#Sweet, not what we have is a tibble with each OTU at all taxonomy levels and the rel-abundance for that for each mouse

#Now, that's a lot of data, but we really only care about one Taxonomy level @ a time

#So, let's filter out Phyla first since it's the largest that gives us useful info
phylum_avg_abund <- otu_rel_abund%>%
  filter(level=="phylum")%>%
  group_by(abx, timepoint, trial,  Group, taxon)%>%  #okay, here we're grouping #this is IMPORTANT: we first need to group by mouse (i.e., Group)
  summarize(rel_abund = sum(rel_abund), .groups = "drop")%>%
  #Now we're grouping to get our average relative abundance, what you want to group by here is whatever major grouping of mice you care about
  group_by(abx, timepoint, trial, taxon)%>%
  summarize(mean_rel_abund = mean(rel_abund), .groups = "drop")

#IDK about you but I like to see % not decimals...
#So just take the abund and * 100 and make a new column for this 
phylum_avg_abund$per_abund<- phylum_avg_abund$mean_rel_abund*100
#tell me what are the unique phylum 
unique(phylum_avg_abund$taxon)
#oops... things look weird with some ""extra"" quotes... let's remove those
phylum_avg_abund <- phylum_avg_abund%>%
  mutate(taxon = gsub("\"", "", taxon))%>%
  mutate(taxon = ifelse(taxon == "Bacteria_unclassified", "Other", taxon))%>%
  mutate(taxon = ifelse(taxon == "Lentisphaerae", "Other", taxon))


#what's the unique ones now
unique(phylum_avg_abund$taxon)

#tell R waht order I want the phyla to be
phylum_avg_abund$taxon <- factor(phylum_avg_abund$taxon, levels = c("Bacteroidetes", "Actinobacteria","Firmicutes", "Verrucomicrobia", "Other", "Proteobacteria"))

#make a plot
phy <- phylum_avg_abund %>%
  ggplot(aes(x = abx, y = per_abund, fill = taxon)) +
  facet_grid(trial~timepoint, scales = "free") + 
  geom_col()+
  scale_fill_manual(values = phy.colors)+
  labs(x = NULL,
       y = "Mean Proportional Relative Abundance (%)") +
  guides(fill = guide_legend(title = "Phylum")) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust= 0.5),
    #customize facet title
    strip.text = element_text(face = "bold", color = "black", size = 12),
    #Customize the background color of the facet titles
    strip.background = element_rect(fill= "lightgrey"),
    # Bold the legend titles
    legend.text = element_text(face = "plain"),
    legend.title = element_text(face = "bold"),
    # Bold the X and Y axis titles
    axis.title.x = element_text(face = "plain", size = 11, vjust = 1),
    axis.title.y = element_text(face = "plain", size = 11, vjust =1),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_line(color = "white"),
    # Adjust legend position
    legend.position = "right", # Place legend at the bottom
    legend.direction = "vertical", # Display legend horizontally
  )

phy#print it 

int <- ggplotly(phy)#make it interactive for funsies

#let's assign some colors for each unique phylum 
phy.colors <- c("Bacteroidetes" = "#486493" , 
                "Actinobacteria" = "gold",
                "Verrucomicrobia" = "#CB89D7", 
                "Proteobacteria" = "#C03F3F", 
                "Other" ="grey", 
                "Firmicutes" = "#619367")




all_abund <- otu_rel_abund%>%
  filter(level %in% c("phylum", "class", "order", "family"))%>%
  group_by(abx, timepoint, trial,  Group, taxon, level)%>%  #okay, here we're grouping #this is IMPORTANT: we first need to group by mouse (i.e., Group)
  summarize(rel_abund = sum(rel_abund), .groups = "drop")%>%
  #Now we're grouping to get our average relative abundance, what you want to group by here is whatever major grouping of mice you care about
  group_by(abx, timepoint, trial, level, taxon)%>%
  summarize(mean_rel_abund = mean(rel_abund), .groups = "drop")

################################################################################
#All Levels 
################################################################################
otu_avg_abund<- otu_rel_abund%>%
  filter(level=="otu")%>%
  group_by(abx, timepoint, trial,  Group, taxon)%>%  #okay, here we're grouping #this is IMPORTANT: we first need to group by mouse (i.e., Group)
  summarize(rel_abund = sum(rel_abund), .groups = "drop")%>%
  #Now we're grouping to get our average relative abundance, what you want to group by here is whatever major grouping of mice you care about
  group_by(abx, timepoint, trial, taxon)%>%
  summarize(mean_rel_abund = mean(rel_abund), .groups = "drop")%>%
  rename(otu = taxon)
#IDK about you but I like to see % not decimals...
#So just take the abund and * 100 and make a new column for this 
otu_avg_abund$per_abund<- otu_avg_abund$mean_rel_abund*100

#make a plot
otu.test <- otu_avg_abund %>%
  ggplot(aes(x = abx, y = per_abund, fill = otu)) +
  facet_grid(trial~timepoint, scales = "free") + 
  geom_col()+
  labs(x = NULL,
       y = "Mean Proportional Relative Abundance (%)") +
  guides(fill = "none") +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust= 0.5),
    #customize facet title
    strip.text = element_text(face = "bold", color = "black", size = 12),
    #Customize the background color of the facet titles
    strip.background = element_rect(fill= "lightgrey"),
    # Bold the legend titles
    legend.text = element_text(face = "plain"),
    legend.title = element_text(face = "bold"),
    # Bold the X and Y axis titles
    axis.title.x = element_text(face = "plain", size = 11, vjust = 1),
    axis.title.y = element_text(face = "plain", size = 11, vjust =1),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_line(color = "white"),
    # Adjust legend position
    legend.position = "top", # Place legend at the bottom
    legend.direction = "horizontal", # Display legend horizontally
  )

otu.test#print it 

#sick! Everything adds up to 100%, so the averages are correctly calculated.

##make reference taxon table##
lookup_tax<- taxonomy%>%
  select(otu, kingdom, phylum,class,order,family,genus)%>%
  mutate(kingdom = gsub("\\s*\\(\\d+\\)$", "", kingdom))%>%
  mutate(phylum = gsub("\\s*\\(\\d+\\)$", "", phylum))%>%
  mutate(class = gsub("\\s*\\(\\d+\\)$", "", class))%>%
  mutate(order = gsub("\\s*\\(\\d+\\)$", "", order))%>%
  mutate(family = gsub("\\s*\\(\\d+\\)$", "", family))%>%
  mutate(genus = gsub("\\s*\\(\\d+\\)$", "", genus))%>%
  mutate(kingdom = gsub("\"", "",kingdom))%>%
  mutate(phylum = gsub("\"", "",phylum))%>%
  mutate(class = gsub("\"", "",class))%>%
  mutate(order = gsub("\"", "",order))%>%
  mutate(family = gsub("\"", "",family))%>%
  mutate(genus = gsub("\"", "",genus))

#now, join the avg abund otu tab with this lookup to give full taxonomy info for each OTU
otu.abund.lab <- inner_join(otu_avg_abund, lookup_tax)%>%
  mutate(phylum = gsub("\"", "", phylum))%>%
  mutate(phylum = ifelse(phylum == "Bacteria_unclassified", "Other", phylum))%>%
  mutate(phylum = ifelse(phylum == "Lentisphaerae", "Other", phylum))
#now tell R what order we want things
unique(otu.abund.lab$phylum)
otu.abund.lab$phylum <- factor(otu.abund.lab$phylum, levels = c("Bacteroidetes", "Actinobacteria","Firmicutes", "Verrucomicrobia", "Other", "Proteobacteria"))
unique(otu.abund.lab$phylum)
unique(otu.abund.lab$abx)
otu.abund.lab <- otu.abund.lab %>%
  mutate(abx = ifelse(abx == "NDC", "No Drug Control (NDC)", abx))%>%
  mutate(abx = ifelse(abx == "IBZ", "Ibezapolstat (IBZ)", abx))%>%
  mutate(abx = ifelse(abx == "FDX", "Fidaxomicin (FDX)", abx))%>%
  mutate(abx = ifelse(abx == "VAN", "Vancomycin (VAN)", abx))%>%
  mutate(abx = ifelse(abx == "MTZ", "Metronidazole (MTZ)", abx))
unique(otu.abund.lab$abx)
otu.abund.lab$abx <- factor(otu.abund.lab$abx, levels = c("No Drug Control (NDC)", "Ibezapolstat (IBZ)", "Fidaxomicin (FDX)", "Vancomycin (VAN)", "Metronidazole (MTZ)"))
unique(otu.abund.lab$timepoint)
otu.abund.lab <- otu.abund.lab %>%
  mutate(timepoint = ifelse(timepoint == "End of Treatment", "ABX Day 10", timepoint))%>%
  mutate(timepoint = ifelse(timepoint == "1 Week Post Change in Diet", "Post Diet ∆", timepoint))%>%
  mutate(timepoint = ifelse(timepoint == "1 Week Post FMT", "Post FMT", timepoint))%>%
  mutate(timepoint = ifelse(timepoint == "2 Days of Treatment", "ABX Day 2", timepoint))
unique(otu.abund.lab$timepoint) 
otu.abund.lab$timepoint <- factor(otu.abund.lab$timepoint, levels = c("Post FMT", "Post Diet ∆", "ABX Day 2", "ABX Day 10"))

phy <- otu.abund.lab %>%
  ggplot(aes(x = timepoint, y = per_abund, fill = phylum)) +
  facet_grid(trial~abx, scales = "free", space = "free_x") + 
  geom_col()+
  scale_fill_manual(values = phy.colors)+
  labs(x = NULL,
       y = "Mean Proportional Relative Abundance (%)", title = "Phylym Level") +
  guides(fill = guide_legend(title = "Phylum")) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust= 0.5),
    #customize facet title
    strip.text = element_text(face = "bold", color = "black", size = 12),
    #Customize the background color of the facet titles
    strip.background = element_rect(fill= "lightgrey"),
    # Bold the legend titles
    legend.text = element_text(face = "plain"),
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", angle=-85, hjust=-0.18),
    axis.title.x = element_text(face = "bold", size = 13, vjust = 1),
    axis.title.y = element_text(face = "plain", size = 11, vjust =1),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_line(color = "white"),
    # Adjust legend position
    legend.position = "right", # Place legend at the bottom
    legend.direction = "vertical", # Display legend horizontally
  )+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), color = "black", fill = NA, inherit.aes = FALSE)


phy

#now paste phylym name to family
unique.fam.per.phy <- aggregate(family ~ phylum, data = otu.abund.lab, FUN = function(x) length(unique(x)))
print(unique.fam.per.phy)
print(sort(unique(otu.abund.lab$fam.name)))
otu.abund.lab$fam.name <- paste("p_",otu.abund.lab$phylum,"_f_",otu.abund.lab$family, sep="")


#now make colors for each unique family (ask ChatGPT)
fam.col <- c("p_Actinobacteria_f_Actinobacteria_unclassified" = "#FFD700",
                   "p_Actinobacteria_f_Bifidobacteriaceae" = "#DAA520",
                   "p_Actinobacteria_f_Coriobacteriaceae" = "#FFED97",
                   "p_Actinobacteria_f_Corynebacteriaceae" = "#B8860B",
                   "p_Actinobacteria_f_Micrococcaceae" = "gold2",
                   "p_Actinobacteria_f_Mycobacteriaceae" = "yellow2",
                   "p_Actinobacteria_f_Nocardiaceae" = "yellow",
                   "p_Actinobacteria_f_Solirubrobacterales_unclassified" = "yellow3",
                   "p_Actinobacteria_f_Streptomycetaceae" = "gold3",
                   "p_Bacteroidetes_f_Bacteroidaceae" = "#486493",
                   "p_Bacteroidetes_f_Bacteroidales_unclassified" = "#1E90FF",
                   "p_Bacteroidetes_f_Porphyromonadaceae" = "#93AFED",
                   "p_Bacteroidetes_f_Prevotellaceae" = "#759AED",
                   "p_Bacteroidetes_f_Rikenellaceae" = "#87CEEB",
                   "p_Firmicutes_f_Acidaminococcaceae" = "#90EE90",
                   "p_Firmicutes_f_Bacillaceae_1" = "#2D5419",
                   "p_Firmicutes_f_Bacillaceae_2" = "#3CB371",
                   "p_Firmicutes_f_Bacillales_unclassified" = "#2E8B57",
                   "p_Firmicutes_f_Bacilli_unclassified" = "#228B22",
                   "p_Firmicutes_f_Clostridia_unclassified" = "#008000",
                   "p_Firmicutes_f_Clostridiaceae_1" = "#32CD32",
                   "p_Firmicutes_f_Clostridiales_unclassified" = "#9ACD32",
                   "p_Firmicutes_f_Enterococcaceae" = "#ADFF2F",
                   "p_Firmicutes_f_Erysipelotrichaceae" = "#7FFF00",
                   "p_Firmicutes_f_Eubacteriaceae" = "#80AF68",
                   "p_Firmicutes_f_Firmicutes_unclassified" = "#A9DEBD",
                   "p_Firmicutes_f_Lachnospiraceae" = "#449A17",
                   "p_Firmicutes_f_Lactobacillaceae" = "#32CD32",
                   "p_Firmicutes_f_Lactobacillales_unclassified" = "#228B22",
                   "p_Firmicutes_f_Leuconostocaceae" = "#008000",
                   "p_Firmicutes_f_Listeriaceae" = "#90EE90",
                   "p_Firmicutes_f_Paenibacillaceae_1" = "forestgreen",
                   "p_Firmicutes_f_Peptococcaceae_1" = "#3CB371",
                   "p_Firmicutes_f_Peptostreptococcaceae" = "#2E8B57",
                   "p_Firmicutes_f_Planococcaceae" = "#228B22",
                   "p_Firmicutes_f_Ruminococcaceae" = "#008000",
                   "p_Firmicutes_f_Streptococcaceae" = "#34714A",
                   "p_Other_f_Bacteria_unclassified" = "grey",
                   "p_Other_f_Victivallaceae" = "grey",
                   "p_Proteobacteria_f_Betaproteobacteria_unclassified" = "#FFA07A",
                   "p_Proteobacteria_f_Bradyrhizobiaceae" = "#FA8072",
                   "p_Proteobacteria_f_Burkholderiales_unclassified" = "#E9967A",
                   "p_Proteobacteria_f_Desulfovibrionaceae" = "#F08080",
                   "p_Proteobacteria_f_Enterobacteriaceae" = "#CD5C5C",
                   "p_Proteobacteria_f_Methylobacteriaceae" = "#FF6347",
                   "p_Proteobacteria_f_Moraxellaceae" = "#FF4500",
                   "p_Proteobacteria_f_Oxalobacteraceae" = "maroon",
                   "p_Proteobacteria_f_Proteobacteria_unclassified" = "#D66243",
                   "p_Proteobacteria_f_Pseudomonadaceae" = "red4",
                   "p_Proteobacteria_f_Sutterellaceae" = "#D0675B",
                   "p_Verrucomicrobia_f_Verrucomicrobiaceae" = "#90362C"
)

fam <- otu.abund.lab %>%
  ggplot(aes(x = timepoint, y = per_abund, fill = fam.name)) +
  facet_grid(trial~abx, scales = "free", space = "free_x") + 
  geom_col()+
  scale_fill_manual(values = fam.col)+
  labs(x = NULL,
       y = "Mean Proportional Relative Abundance (%)", title = "Family Level") +
  guides(fill = guide_legend(title = "Family")) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust= 0.5),
    #customize facet title
    strip.text = element_text(face = "bold", color = "black", size = 12),
    #Customize the background color of the facet titles
    strip.background = element_rect(fill= "lightgrey"),
    # Bold the legend titles
    legend.text = element_text(face = "plain"),
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", angle=-85, hjust=-0.18),
    axis.title.x = element_text(face = "bold", size = 13, vjust = 1),
    axis.title.y = element_text(face = "plain", size = 11, vjust =1),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_line(color = "white"),
    # Adjust legend position
    legend.position = "right", # Place legend at the bottom
    legend.direction = "vertical", # Display legend horizontally
  )+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), color = "black", fill = NA, inherit.aes = FALSE)


fam

#now paste phylym name to order
unique.ord.per.phy <- aggregate(order ~ phylum, data = otu.abund.lab, FUN = function(x) length(unique(x)))
print(unique.ord.per.phy)
print(sort(unique(otu.abund.lab$ord.name)))
otu.abund.lab$ord.name <- paste("p_",otu.abund.lab$phylum,"_o_",otu.abund.lab$order, sep="")

ord.col <- c(
  "p_Actinobacteria_o_Actinobacteria_unclassified" = "#FFD700",
  "p_Actinobacteria_o_Actinomycetales" = "#DADF92",
  "p_Actinobacteria_o_Bifidobacteriales" = "#9DA42F",
  "p_Actinobacteria_o_Coriobacteriales" = "#8C884B",
  "p_Actinobacteria_o_Solirubrobacterales" = "#DAD368",
  "p_Bacteroidetes_o_Bacteroidales" = "#4D6082",
  "p_Firmicutes_o_Bacillales" = "#77A179",
  "p_Firmicutes_o_Bacilli_unclassified" = "#7C986E",
  "p_Firmicutes_o_Clostridia_unclassified" = "#B7D69D",
  "p_Firmicutes_o_Clostridiales" = "#176534",
  "p_Firmicutes_o_Erysipelotrichales" = "#8CCB49",
  "p_Firmicutes_o_Firmicutes_unclassified" = "#D7EEBF",
  "p_Firmicutes_o_Lactobacillales" = "#1EB864",
  "p_Firmicutes_o_Selenomonadales" = "#264D38",
  "p_Other_o_Bacteria_unclassified" = "#808080",
  "p_Other_o_Victivallales" = "#808080",
  "p_Proteobacteria_o_Betaproteobacteria_unclassified" = "#FFA07A",
  "p_Proteobacteria_o_Burkholderiales" = "#FF7F7F",
  "p_Proteobacteria_o_Desulfovibrionales" = "#FF4500",
  "p_Proteobacteria_o_Enterobacteriales" = "#D0675B",
  "p_Proteobacteria_o_Proteobacteria_unclassified" = "#EDBABA",
  "p_Proteobacteria_o_Pseudomonadales" = "#6E2626",
  "p_Proteobacteria_o_Rhizobiales" = "#C2593A",
  "p_Verrucomicrobia_o_Verrucomicrobiales" = "#E5B3FF"
)

ord<- otu.abund.lab %>%
  ggplot(aes(x = timepoint, y = per_abund, fill = ord.name)) +
  facet_grid(trial~abx, scales = "free", space = "free_x") + 
  geom_col()+
  scale_fill_manual(values = ord.col)+
  labs(x = NULL,
       y = "Mean Proportional Relative Abundance (%)", title = "Order Level") +
  guides(fill = "none") +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust= 0.5),
    #customize facet title
    strip.text = element_text(face = "bold", color = "black", size = 12),
    #Customize the background color of the facet titles
    strip.background = element_rect(fill= "lightgrey"),
    # Bold the legend titles
    legend.text = element_text(face = "plain"),
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", angle=-85, hjust=-0.18),
    axis.title.x = element_text(face = "bold", size = 13, vjust = 1),
    axis.title.y = element_text(face = "plain", size = 11, vjust =1),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_line(color = "white"),
    # Adjust legend position
    legend.position = "right", # Place legend at the bottom
    legend.direction = "vertical", # Display legend horizontally
  )+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), color = "black", fill = NA, inherit.aes = FALSE)


ord

#now paste phylym name to class
unique.cls.per.phy <- aggregate(class ~ phylum, data = otu.abund.lab, FUN = function(x) length(unique(x)))
print(unique.cls.per.phy)
print(sort(unique(otu.abund.lab$cls.name)))
otu.abund.lab$cls.name <- paste("p_",otu.abund.lab$phylum,"_o_",otu.abund.lab$class, sep="")

cls.col <- c(
  "p_Actinobacteria_o_Actinobacteria" = "gold",
  "p_Bacteroidetes_o_Bacteroidia" = "#4D6082",
  "p_Firmicutes_o_Bacilli" = "#77A179",
  "p_Firmicutes_o_Clostridia" = "#176534",
  "p_Firmicutes_o_Erysipelotrichia" = "#8CCB49",
  "p_Firmicutes_o_Firmicutes_unclassified" = "#D7EEBF",
  "p_Firmicutes_o_Negativicutes" = "#619367",
  "p_Other_o_Bacteria_unclassified" = "#808080",
  "p_Other_o_Lentisphaeria" = "#808080",
  "p_Proteobacteria_o_Alphaproteobacteria" = "#C2593A",
  "p_Proteobacteria_o_Betaproteobacteria" = "#FF7F7F",
  "p_Proteobacteria_o_Deltaproteobacteria" = "#B63E2B",
  "p_Proteobacteria_o_Gammaproteobacteria" = "#8A1908",
  "p_Proteobacteria_o_Proteobacteria_unclassified" = "#F99485",
  "p_Verrucomicrobia_o_Verrucomicrobiae" = "#E5B3FF"
)

cls<- otu.abund.lab %>%
  ggplot(aes(x = timepoint, y = per_abund, fill = cls.name)) +
  facet_grid(trial~abx, scales = "free", space = "free_x") + 
  geom_col()+
  scale_fill_manual(values = cls.col)+
  labs(x = NULL,
       y = "Mean Proportional Relative Abundance (%)", title = "Class Level") +
  guides(fill = guide_legend(title = "Class")) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust= 0.5),
    #customize facet title
    strip.text = element_text(face = "bold", color = "black", size = 12),
    #Customize the background color of the facet titles
    strip.background = element_rect(fill= "lightgrey"),
    # Bold the legend titles
    legend.text = element_text(face = "plain"),
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", angle=-85, hjust=-0.18),
    axis.title.x = element_text(face = "bold", size = 13, vjust = 1),
    axis.title.y = element_text(face = "plain", size = 11, vjust =1),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_line(color = "white"),
    # Adjust legend position
    legend.position = "right", # Place legend at the bottom
    legend.direction = "vertical", # Display legend horizontally
  )+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), color = "black", fill = NA, inherit.aes = FALSE)


cls



