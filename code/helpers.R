

# convert a list of dataframe into a long dataframe
list_to_df <- function(mylist){
  
  # make a vector of row ids that correspond to the list names
  rowid.indx <- lapply(mylist, function(x) dim(x)[1])
  sourceVec.list <- list()
  for(i in 1:length(rowid.indx)){
    sourceName <- names(rowid.indx)[i]
    numRows <- rowid.indx[[i]]
    sourceVec.list[[i]] <- rep(sourceName, numRows)
  }
  rowVec <- unlist(sourceVec.list)
  
  # combine into df
  df <- data.frame(do.call(rbind, mylist))
  df$source <- rowVec
  
  return(df)
}

# source a bunch of files from the same folder; code taken from R help -> source()
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

# create a new folder
add.dir <- function(dir){
  if(!dir.exists(dir)) dir.create(dir)
}

# extend x/y lims to make the MDS plot a square
square_lims <- function(df){
  
  df %>%
    filter(!is.na(MDS1)) %>%
    filter(!is.na(MDS2)) -> df
  
  curr.xlims <- range(df$MDS1)
  curr.ylims <- range(df$MDS2)
  xlen <- range(df$MDS1)[2] - range(df$MDS1)[1]
  ylen <- range(df$MDS2)[2] - range(df$MDS2)[1]
  max.len <- max(c(xlen, ylen))
  if(max.len == xlen){
    addthis <- max.len - ylen
    curr.ylims[2]<- curr.ylims[2] + addthis # extend ylen
  }else{
    addthis <- max.len - xlen
    curr.xlims[2]<- curr.xlims[2] + addthis # extend xlen
  }
  
  result <- list(xlims = curr.xlims,
                 ylims = curr.ylims)
  return(result)
  
}

# custom pca biplot
ggbiplot_prcomp <- function(pca, orig.df, select.pcs){
  
  # make df for the site scores
  scrs <- data.frame(orig.df, 
                     PC1 = pca$x[,select.pcs[1]], 
                     PC2 = pca$x[,select.pcs[2]])
  
  # make df for the loadings (vectors)
  spp.scrs<- data.frame(var = row.names(pca$rotation),
                        PC1 = pca$rotation[,select.pcs[1]], 
                        PC2 = pca$rotation[,select.pcs[2]])
  
  # make axis labels
  pca.summ <- summary(pca)
  pc1.import <- round(pca.summ$importance[2,select.pcs[1]], 
                      digits = 4) * 100
  pc2.import <- round(pca.summ$importance[2,select.pcs[2]], 
                      digits = 4) * 100
  xlablabel <- paste0("PC", select.pcs[1]," ",pc1.import,"%")
  ylablabel <- paste0("PC", select.pcs[2]," ",pc2.import,"%")
  
  # plot
  mult <- 10 #multiplier for the arrows and text for envfit
  mult.text <- 10.4
  
  p <- ggplot(scrs) +
    geom_point(mapping = aes(x = PC1, y = PC2)) +
    coord_fixed() + ## need aspect ratio of 1!
    xlab(xlablabel) + 
    ylab(ylablabel)
  p <- p + geom_segment(data = spp.scrs, color = "red",
                        aes(x = 0, xend = mult*PC1, 
                            y = 0, yend = mult*PC2),
                        arrow = arrow(length = unit(0.2, "cm"))) +
    geom_text(data = spp.scrs, color = "red",
              aes(x = mult.text*PC1, 
                                   y = mult.text*PC2, 
                                   label = var), 
              size = 4) 
  
  return(p)
  
}
