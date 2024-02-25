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
  library(ggforce)
}


################################################################################
#Using mothur output files to make nmds plots & beta diversity analysis        
################################################################################

#First, we need to read in the mothur output file and format it to work with us in R
#you will read in your "final.opti_mcc.shared" file from mothur
count.otu<- read_tsv("/Users/trentonwolfe/UofH-Collab/IBZ Manuscript/mothur/mothur.out/final.opti_mcc.shared")%>%
  select(Group, starts_with("Otu"))%>%
  pivot_longer(-Group)
count.otu#this will output a tibble that shows the sample ID, an Otu, and the count that that specific OTU appeared within the seqs from that sample
#we now will read in the metadata2 that we created in the set-up script and join that with this tibble
lab.count.otu <- inner_join(count.otu,metadata)#this should automatically join by 'Group', which is the sampleID

#now, lets make it a table by mutating and pivoting the above df to make each unique otu its own column witht the data in the rows the count for each otu for each sample ('Group')
otu.tbl<- lab.count.otu%>%
  group_by(Group)%>%
  mutate(N = sum(value))%>%
  ungroup()%>%
  pivot_wider(names_from = "name", values_from = "value", values_fill = 0)

#now the next two steps, we need to make two df that seperate what is metadata and what is data to make a distance matrix

#first, let's pull out all things that aren't metadata
meta.tbl <- otu.tbl%>%
  select(Group, mouse, timepoint, day, `abx-grp`, treatment, abx, sex, trial, food, N, ellipse)
#okay, now we can basically subtract ALL of these (plus the label thing mothur adds) and pipeline in avgdist to make a distance matrix
dist.mat <- otu.tbl%>%
  select(-mouse, -timepoint, -day, -`abx-grp`, -treatment, -abx, -sex, -trial, -food, -N, -ellipse)%>%
  column_to_rownames("Group")%>%
  avgdist(sample=1018)#I set this based on what Pat Schloss used. It does kick out a few samples, but keeps the good ones. Theres 5 samples that were below this #, which is fine.
?metaMDS
#coolio, now let's set the seed and then make an nmds tibble from this distance matrix
set.seed(1)
nmds.tbl <- metaMDS(dist.mat)%>%
  scores()%>%
  as_tibble(rownames="Group")

#now, let's rejoin in our metadata variables we pulled earlier and rejoin them 
meta.nmds <- inner_join(meta.tbl, nmds.tbl)

#AMAZING: you just now made nmds coordinates to make an nmds plot with the data and now have everything all in one df to make it very easy to plot and manipulate with ggplot
#first, let's tell R what order things happened to make things nice and easy for the casual reader to understand the data

unique(meta.nmds$timepoint)
meta.nmds$timepoint <- factor(meta.nmds$timepoint, levels = c("1 Week Post FMT", "1 Week Post Change in Diet", "2 Days of Treatment", "End of Treatment"))


#time to start plotting
nmds_plot<- meta.nmds%>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=abx, shape=trial, group=ellipse))+
  geom_point(size=5)+
  facet_grid(timepoint ~ ., scales = "fixed")+
  scale_color_manual(values= abx.col)+
  guides(color= guide_legend(title = "Treatment"))+
  guides(shape= guide_legend(title = "Trial"))+
  labs(title = "",
       x = "NMDS1", y = "NMDS2")+theme(text = element_text(family = "Arial"))+
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
  )+
  scale_x_continuous(breaks = seq(-0.7, 0.4, by = 0.1), limits = c(-0.7, 0.4),
                     labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(breaks = seq(-0.6, 0.6, by = 0.1), limits = c(-0.6, 0.6),
                     labels = scales::number_format(accuracy = 0.1)) +
  geom_mark_ellipse()+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), color = "black", fill = NA, inherit.aes = FALSE)





#print the plot by simply providing its name
nmds_plot   

#now let's add some ellipse around each group with a 95% CI
nmds_treatellipse<- nmds_plot+
  stat_ellipse(aes(group=ellipse), type= "t", level = 0.9)

nmds_treatellipse
#for this specific dataset, we didn't have large arms for each ABX in each trial... thus the ellipses weren't generated with too small n
#however, because we're fancy here, let's animate the plot to see the temporal microbiome changes 
#time to start plotting :) (basically, the same as above, but without the facet_wrap/facet_grid)
nmds_plot_int<- meta.nmds%>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=abx, shape=trial, group=abx))+
  geom_point(size=5)+
  scale_color_manual(values= abx.col)+
  guides(color= guide_legend(title = "Group"))+
  guides(shape= guide_legend(title = "Trial"))+
  labs(title = "",
       x = "NMDS1", y = "NMDS2")+theme(text = element_text(family = "Arial"))+
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust= 0.5),
    #customize facet title
    strip.text = element_text(face = "bold", color = "black", size = 12),
    #Customize the background color of the facet titles
    strip.background = element_rect(fill= "white"),
    # Bold the legend titles
    legend.text = element_text(face = "plain"),
    legend.title = element_text(face = "bold"),
    # Bold the X and Y axis titles
    axis.title.x = element_text(face = "bold", size = 12, vjust = 1),
    axis.title.y = element_text(face = "bold", size = 12),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "lightgrey"),
    panel.grid.minor = element_line(color = "white"),
    # Adjust legend position
    legend.position = "top", # Place legend at the bottom
    legend.direction = "horizontal" # Display legend horizontally
  )+coord_cartesian(ylim = ylim, xlim= xlim)+ geom_mark_ellipse()

?geom_mark_ellipse
#now, to animate, this is a very straightforward one function... thanks R gods... 
gif <- nmds_plot_int+transition_states(timepoint, transition_length = 1000)+
  ease_aes('linear')+
  labs(subtitle = "Timepoint: {closest_state}")

#to view, we simply just do:
gif


#Neat, huh? Now, let's write a file to save this animation for later use 
anim_save("nmds.animation.gif", animation = gif)
#the "" will be the file name saved in the current working directory of where your R project is stored.

#that's all folks... 







