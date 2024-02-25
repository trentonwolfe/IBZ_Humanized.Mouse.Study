#PERMANOVA ANALYSIS FOR NMDS PLOTS
#To do some PERMANOVA analysis for the NMDS plots, we first need a df that has just OTU counts for each OTU and the rownames as the Group
otu.stat <- otu.tbl%>%
  select(-mouse, -timepoint, -day, -`abx-grp`, -treatment, -abx, -sex, -trial, -food, -ellipse, -N)%>%
  column_to_rownames(var = "Group")
#NEAT!




#So now let's make a loop that iterates through each day
for(i in unique(metadata$timepoint)){
  #Make a sub df for each i 
  sub<- metadata$timepoint == i
  #Print a line to make it easy to differentiate between each output
  print("_______________________________________")
  #Tell me what day this analysis is for
  print(paste("At [", i, "] PERMANOVA Test of Trial (FMT Donor) Explanation of Variation"))
  #Analyze the OTU count df to ask if they are different on an NMDS plot because of what vendor they originated from (I have this in the metadata file in the experiment column)
  res2<- adonis2(otu.stat[sub,]~metadata$abx[sub], method = "bray")
  #Print the PERMANOVA results
  print(res2)
  #Close out with a line 
  print("_______________________________________")
}


##permanova code on NMDS plot (Fig 2 F)
library(vegan)
o <- order(metadata$day)

otu <- otu.stat[o,]
map <- map[o,]

map$day

all(row.names(otu.stat) == meta.tbl$Group)

mice <-unique(map$`Mouse-ID`)
days <- unique(map$day)

res <- matrix(ncol = 4, nrow = 32)


for(i in 1 : length(mice)){
  
  otu.i <- otu.stat[meta.tbl$Group == mice[i], ]
  map.i <- map[map$`Mouse-ID` == mice[i], ]
  
  map.i$day
  
  aa <- as.matrix(vegdist(otu.i))[1,]
  
  res[i,   days %in% map.i$day] <- aa
}

##
plot(res[1, ], type = "l", ylim = range(res, na.rm = T))
for(i in 2 : length(mice)){
  
  points(res[i, ], type = "l")
}


##
library(vegan)
#significantly different for each day by means of treatment given during ABX? 
for(i in unique(meta.tbl$day)){
  
  p <- meta.tbl$day == i
  print(paste("For Day:", i, "of the Experiment"))
  res <- adonis2(otu.stat[p,] ~ meta.tbl$abx[p])
  print(res)
}

#significantly different for each day by menas of trial included in? 
for(i in unique(meta.tbl$day)){
  
  p <- meta.tbl$day == i
  print("-------------------------")
  print(paste("For Day:", i, "of the Experiment"))
  res <- adonis2(otu.stat[p,] ~ meta.tbl$trial[p])
  print(res)
  print("--------End of Test--------")
  
}




