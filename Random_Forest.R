#Random Forest Analysis. Primarily Written by Nick V. Pinkham

remove_rare <- function( table , cutoff_pro ) {
  
  table <- t(table)# transpose to keep "tidy" ; easier that rewriting function...
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) )
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  res <- t(table [ row2keep , , drop=F ])
  res <- as.data.frame(res)
  
  return(res )
}

map$Discription <- map$antibiotic.group


antibiotics <- unique(map$Discription)

taxonomy.pick <- taxonomy[match(colnames(otu),taxonomy$otu) , ]
taxonomy.pick$family <- gsub("\\(\\d+\\)", "", taxonomy.pick$family)

all(colnames(otu) == taxonomy.pick$otu)
#ONLY RUN ONCE!
colnames(otu) <- paste(colnames(otu), taxonomy.pick$family)


#make a function that does RF analysis and prints the plot
pairwiseRF <- function(mat, map, dis_1 = "T 0", dis_2 = "B 7",
                       num_of_col = 50, Accuray.cutoff = 0,
                       Factor = 1, plot = T){
  # report mats that are different between two treatments based on RF analysis
  set.seed(42)
  mat.pick <- mat[map$Discription %in% c(dis_1, dis_2 ), ]
  map.pick <- map[map$Discription %in% c(dis_1, dis_2 ), ]
  
  #mat.pick$Group <- NULL
  
  mat.pick.rare_removed <- remove_rare(table = mat.pick, cutoff_pro = 0.2)
  
  spliter <- as.factor(as.character(map.pick$Discription))
  
  fit <- randomForest::randomForest(x = mat.pick.rare_removed,
                                    y = spliter,
                                    importance=TRUE,
                                    proximities=TRUE,
                                    ntree = 5001)
  RF.sig <- rfUtilities::rf.significance(x = fit,
                                         xdata = mat.pick.rare_removed,
                                         nperm = 101,
                                         ntree = 501)
  print(fit)
  print(RF.sig)
  
  if(RF.sig$pValue < 0.05){
    
    fit_imp <- as.data.frame( fit$importance )
    fit_imp$features <- rownames(fit_imp )
    fit_imp_sorted <- dplyr::arrange( fit_imp  , desc(fit_imp$MeanDecreaseGini)  )
    colors <- vector(length = ncol(mat.pick.rare_removed))
    
    for(j in 1:ncol(mat.pick.rare_removed)){
      i <- fit_imp_sorted$features[j]
      t1.mean <- mean(mat.pick[map.pick$Discription == dis_1, which(colnames(mat.pick) == i)])
      t2.mean <- mean(mat.pick[map.pick$Discription == dis_2, which(colnames(mat.pick) == i)])
      if( t1.mean >  t2.mean){
        colors[j] <- "cadetblue3"
      }else{
        colors[j] <- "orange3"
      }
    }
    
    if(plot){
      par(mar = c(12, 4, 2, 2) + 0.1)
      barplot(fit_imp_sorted$MeanDecreaseGini[1:num_of_col],
              names.arg = fit_imp_sorted[1:num_of_col,"features"],
              ylab= "OTU Importance",
              las = 2,
              col = colors,
              cex.names = 0.8,
              main= paste("RF Comparison of", dis_1,"vs", dis_2))
      
      legend("topright",
             legend = paste("Elevated in:", c(dis_1, dis_2)),
             fill = c("cadetblue3", "orange3"),
             border = "black",
             bty = "n",
             xjust=0.5,
             yjust=0.5)
      
    }
    
    res <- fit_imp_sorted$MeanDecreaseGini
    names(res) <- fit_imp_sorted$features
    
    if(plot){
      abline(h = median(res) + sd(res) * Factor,
             col = 2, lty = 2)
      
    }
    res.pick <- res[res > median(res) + (sd(res) * Factor)]
    colors.pick <- colors[res > median(res) + (sd(res) * Factor)]
    
    mat.bloom <- res.pick[colors.pick == "orange3"]
    mat.bust <- res.pick[colors.pick == "cadetblue3"]
    
    Result <- list(  mat.bust, mat.bloom)
    names(Result) <- c(dis_1, dis_2)
    
    return(Result)
  }else{
    return("VARIABLE p val > 0.05")
  }
}





#loop thorugh each pariwise comparison and run the above function 
for(i in 1 : length(antibiotics)){
  for(j in 1 : i){
    if(i != j){
      res <- pairwiseRF(mat = otu, 
                        map, 
                        dis_1 = antibiotics[i], 
                        dis_2 = antibiotics[j])
      
      res[[1]]
      res[[2]]
      names(res)
      
      names(res$Fidaxomicin)
      
    }
  }
}


#example of how to code and run a desired comparison with function made above
pairwiseRF(mat = otu, map, dis_1 = "Ibezapolstat", dis_2 = "ND Control")


