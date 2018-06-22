plot_ordered_bar2<-function (physeq, x = "Sample", 
                             y = "Abundance", 
                             fill = NULL, 
                             leg_size = 0.5,
                             title = NULL) {
  require(ggplot2)
  require(phyloseq)
  require(plyr)
  require(grid)
  bb <- psmelt(physeq)
  
  samp_names <- aggregate(bb$Abundance, by=list(bb$Sample), FUN=sum)[,1]
  .e <- environment()
  bb[,fill]<- factor(bb[,fill], rev(sort(unique(bb[,fill])))) #fill to genus
  
  
  bb<- bb[order(bb[,fill]),] # genus to fill
  p = ggplot(bb, aes_string(x = x, y = y, 
                            fill = fill), 
             environment = .e, ordered = FALSE)
  
  
  p = p +geom_bar(stat = "identity", 
                  position = "stack") 
  
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  
  p = p + guides(fill = guide_legend(override.aes = list(colour = NULL), reverse=TRUE)) + 
    theme(legend.key = element_rect(colour = "black")) 
  
  p = p + theme(legend.key.size = unit(leg_size, "cm"))
  
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}