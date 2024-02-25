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
  library(multcompView)
  
}


################################################################################
#Using metadata2 df that was made in the set-up script for diversity metric analysis      
################################################################################
metadata2$abx <- factor(metadata2$abx, levels = c("NDC", "IBZ", "FDX", "VAN", "MTZ"))

meta.nmds$timepoint <- factor(meta.nmds$timepoint, levels = c("1 Week Post FMT", "1 Week Post Change in Diet", "2 Days of Treatment", "End of Treatment"))
#Simpson's Index
for(i in unique(metadata2$timepoint)) {
  sub.i <- metadata2[metadata2$timepoint == i, ]  # Subset data for the current time point

  # Simpson's Index plot
  simpson.plt <- ggplot(sub.i, aes(x=abx, y=simpson, fill = abx, color=abx))+
    geom_boxplot(alpha=0.6)+
    scale_fill_manual(values= abx.col)+
    scale_color_manual(values = abx.col) +
    guides(fill="none")+guides(color="none")+
    labs(title = i, 
         y = "Simpson's Index", x = "") +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(face = "bold", size = 14, hjust= 0.5),
      strip.text = element_text(face = "bold", color = "black", size = 12),
      strip.background = element_rect(fill= "lightgrey"),
      legend.text = element_text(face = "plain"),
      legend.title = element_text(face = "bold"),
      axis.title.x = element_text(face = "bold", size = 11, vjust = 1),
      axis.title.y = element_text(face = "bold", size = 11),
      axis.text.x = element_text(face="bold", size = 11, angle = 18,vjust = 1, hjust = 1),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "lightgrey"),
      panel.grid.minor = element_line(color = "white"),
      legend.position = "right",
      legend.direction = "vertical"
    )+ scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1))+
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), color = "black", fill = NA, inherit.aes = FALSE)
  
  print(paste("Timepoint", i, "of Experiment"))
  print(summary(aov(simpson ~ abx, data = sub.i)))  # Analysis of variance
  
  
  
  # Tukey's Test
  tukey <- TukeyHSD(aov(simpson ~ abx, data = sub.i))
  print(tukey)
  
  # Compact letter display
  cld <- multcompLetters4(aov(simpson ~ abx, data = sub.i), tukey)
  print(cld)
  
  # Table with factors and 3rd quantile
  Tk <- sub.i %>%
    group_by(abx) %>%
    summarise(mean = mean(simpson), quant = quantile(simpson, probs = 0.75)) %>%
    arrange(desc(mean))
  cld <- as.data.frame.list(cld$abx)
  Tk$cld <- cld$Letters
  
  p2 <- simpson.plt + 
    geom_text(data = Tk, aes(label = cld, x = abx, y = 1), color = "black", size = 8)
  
  print(p2)
}

#Shannon's Index
for(i in unique(metadata2$timepoint)) {
  sub.i <- metadata2[metadata2$timepoint == i, ]  # Subset data for the current time point
  
  # Simpson's Index plot
  shannon.plt <- ggplot(sub.i, aes(x=abx, y=shannon, fill = abx, color=abx))+
    geom_boxplot(alpha=0.6)+
    scale_fill_manual(values= abx.col)+
    scale_color_manual(values = abx.col) +
    guides(fill="none")+guides(color="none")+
    labs(title = i, 
         y = "Shannon's Index", x = "") +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(face = "bold", size = 14, hjust= 0.5),
      strip.text = element_text(face = "bold", color = "black", size = 12),
      strip.background = element_rect(fill= "lightgrey"),
      legend.text = element_text(face = "plain"),
      legend.title = element_text(face = "bold"),
      axis.title.x = element_text(face = "bold", size = 11, vjust = 1),
      axis.title.y = element_text(face = "bold", size = 11),
      axis.text.x = element_text(face="bold", size = 11, angle = 18,vjust = 1, hjust = 1),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "lightgrey"),
      panel.grid.minor = element_line(color = "white"),
      legend.position = "right",
      legend.direction = "vertical"
    )+scale_y_continuous(breaks = seq(0, 4, by = 0.5), limits = c(0, 4))+
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), color = "black", fill = NA, inherit.aes = FALSE)
  
  print(paste("Timepoint", i, "of Experiment"))
  print(summary(aov(shannon ~ abx, data = sub.i)))  # Analysis of variance
  
  
  
  # Tukey's Test
  tukey <- TukeyHSD(aov(shannon ~ abx, data = sub.i))
  print(tukey)
  
  # Compact letter display
  cld <- multcompLetters4(aov(shannon ~ abx, data = sub.i), tukey)
  print(cld)
  
  # Table with factors and 3rd quantile
  Tk <- sub.i %>%
    group_by(abx) %>%
    summarise(mean = mean(shannon), quant = quantile(shannon, probs = 0.75)) %>%
    arrange(desc(mean))
  cld <- as.data.frame.list(cld$abx)
  Tk$cld <- cld$Letters
  
  p2 <- shannon.plt + 
    geom_text(data = Tk, aes(label = cld, x = abx, y = 4), color = "black", size = 8)
  
  print(p2)
}
metadata2$sobs...5
#Sobs (Richness) Index
for(i in unique(metadata2$timepoint)) {
  sub.i <- metadata2[metadata2$timepoint == i, ]  # Subset data for the current time point
  
  # Simpson's Index plot
  sobs.plt <- ggplot(sub.i, aes(x=abx, y=sobs...5, fill = abx, color=abx))+
    geom_boxplot(alpha=0.6)+
    scale_fill_manual(values= abx.col)+
    scale_color_manual(values = abx.col) +
    guides(fill="none")+guides(color="none")+
    labs(title = i, 
         y = "Observed Species Richness (Sobs)", x = "") +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(face = "bold", size = 14, hjust= 0.5),
      strip.text = element_text(face = "bold", color = "black", size = 12),
      strip.background = element_rect(fill= "lightgrey"),
      legend.text = element_text(face = "plain"),
      legend.title = element_text(face = "bold"),
      axis.title.x = element_text(face = "bold", size = 11, vjust = 1),
      axis.title.y = element_text(face = "bold", size = 11),
      axis.text.x = element_text(face="bold", size = 11, angle = 18,vjust = 1, hjust = 1),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "lightgrey"),
      panel.grid.minor = element_line(color = "white"),
      legend.position = "right",
      legend.direction = "vertical"
    )+scale_y_continuous(breaks = seq(0, 160, by = 10), limits = c(0, 160))+
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), color = "black", fill = NA, inherit.aes = FALSE)
  
  print(paste("Timepoint", i, "of Experiment"))
  print(summary(aov(sobs...5 ~ abx, data = sub.i)))  # Analysis of variance
  
  
  
  # Tukey's Test
  tukey <- TukeyHSD(aov(sobs...5 ~ abx, data = sub.i))
  print(tukey)
  
  # Compact letter display
  cld <- multcompLetters4(aov(sobs...5 ~ abx, data = sub.i), tukey)
  print(cld)
  
  # Table with factors and 3rd quantile
  Tk <- sub.i %>%
    group_by(abx) %>%
    summarise(mean = mean(sobs...5), quant = quantile(sobs...5, probs = 0.75)) %>%
    arrange(desc(mean))
  cld <- as.data.frame.list(cld$abx)
  Tk$cld <- cld$Letters
  
  p2 <- sobs.plt + 
    geom_text(data = Tk, aes(label = cld, x = abx, y = 160), color = "black", size = 8)
  
  print(p2)
}

##beta dispersion analysis 
bd<- betadisper(dist.mat, meta.tbl$abx)
anova(bd)
permutest(bd)



#make a df that just have the distances from betadisper function from above
bd.df<- as.data.frame(bd$distances)
#make it so I can innerjoin metadata to the df to make it a usable df for gglplot
bd.df<- bd.df%>%
  rownames_to_column(var="Group")%>%
  rename(betadisper = 'bd$distances')
bd.df<- inner_join(bd.df,metadata2)

#Tell R what order I want
bd.df$abx <- factor(bd.df$abx, levels = c("NDC", "IBZ", "FDX", "VAN", "MTZ"))

bd.df$timepoint <- factor(bd.df$timepoint, levels = c("1 Week Post FMT", "1 Week Post Change in Diet", "2 Days of Treatment", "End of Treatment"))

#Betadispersion
for(i in unique(bd.df$timepoint)) {
  sub.i <- bd.df[bd.df$timepoint == i, ]  # Subset data for the current time point
  
  # Simpson's Index plot
  beta.plt <- ggplot(sub.i, aes(x=abx, y=betadisper, fill = abx, color=abx))+
    geom_boxplot(alpha=0.6)+
    scale_fill_manual(values= abx.col)+
    scale_color_manual(values = abx.col) +
    guides(fill="none")+guides(color="none")+
    labs(title = i, 
         y = "Beta Dispersion", x = "") +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(face = "bold", size = 14, hjust= 0.5),
      strip.text = element_text(face = "bold", color = "black", size = 12),
      strip.background = element_rect(fill= "lightgrey"),
      legend.text = element_text(face = "plain"),
      legend.title = element_text(face = "bold"),
      axis.title.x = element_text(face = "bold", size = 11, vjust = 1),
      axis.title.y = element_text(face = "bold", size = 11),
      axis.text.x = element_text(face="bold", size = 11, angle = 18,vjust = 1, hjust = 1),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "lightgrey"),
      panel.grid.minor = element_line(color = "white"),
      legend.position = "right",
      legend.direction = "vertical"
    )+scale_y_continuous(breaks = seq(0, 0.8, by = 0.1), limits = c(0, 0.8))+
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), color = "black", fill = NA, inherit.aes = FALSE)
  
  print(paste("Timepoint", i, "of Experiment"))
  print(summary(aov(betadisper ~ abx, data = sub.i)))  # Analysis of variance
  
  
  
  # Tukey's Test
  tukey <- TukeyHSD(aov(betadisper ~ abx, data = sub.i))
  print(tukey)
  
  # Compact letter display
  cld <- multcompLetters4(aov(betadisper ~ abx, data = sub.i), tukey)
  print(cld)
  
  # Table with factors and 3rd quantile
  Tk <- sub.i %>%
    group_by(abx) %>%
    summarise(mean = mean(betadisper), quant = quantile(betadisper, probs = 0.75)) %>%
    arrange(desc(mean))
  cld <- as.data.frame.list(cld$abx)
  Tk$cld <- cld$Letters
  
  p2 <- beta.plt + 
    geom_text(data = Tk, aes(label = cld, x = abx, y = 0.8), color = "black", size = 8)
  
  print(p2)
}






