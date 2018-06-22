#Phyloseq Stuff
detach_package <- function(pkg, character.only = FALSE)
   {
         if(!character.only)
           {
                 pkg <- deparse(substitute(pkg))
             }
         search_item <- paste("package", pkg, sep = ":")
         while(search_item %in% search())
           {
                 detach(search_item, unload = TRUE, character.only = TRUE)
             }
     }
detach_package(ggtree)
library("phyloseq")
library("ggplot2")
library("plyr")

setwd('C:/Users//Tim Hollibaugh/Google Drive/Hollibaugh Lab/ML16S/qiimeML/Denoised/or_us61_silva/')

#this only imports 238 OTUs (since 5 could not be put into tree)
ML16S_final_biom_w_tree_w_metadata_w_refseq<-import_biom('otu_table_mc2_w_tax_noChloroplast_2012only_noSingletons_updated_taxonomy_w_metadata2.biom', treefilename='FINAL_biom_filtered_rep_set.tre', refseqfilename='FINAL_biom_filtered_rep_set_edited.fna', parseFunction=parse_taxonomy_greengenes)

#this imports 243 OTUs
ML16S_final_biom_w_metadata_w_refseq<-import_biom('otu_table_mc2_w_tax_noChloroplast_2012only_noSingletons_updated_taxonomy_w_metadata2.biom', refseqfilename='FINAL_biom_filtered_rep_set_edited.fna', parseFunction=parse_taxonomy_greengenes)

#Distance
theme_set(theme_bw())
dist_methods <- unlist(distance("list"))
print(dist_methods)
dist_methods[(1:3)]
dist_methods <- dist_methods[-(1:3)]
dist_methods = dist_methods[-which(dist_methods == "ANY")]
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for (i in dist_methods) {
  # Calculate distance matrix
  iDist <- distance(ML16S_final_biom_w_metadata_w_refseq, method = i)
  # Calculate ordination
  iMDS <- ordinate(ML16S_final_biom_w_metadata_w_refseq, "MDS", distance = iDist)
  ## Make plot Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(ML16S_final_biom_w_metadata_w_refseq, iMDS, color = "Temperature", shape = "Depth")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep = ""))
  # Save the graphic to file.
  plist[[i]] = p
}
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color = Temperature, shape = Depth))
p = p + geom_point(size = 3, alpha = 0.5)
p = p + facet_wrap(~distance, scales = "free")
p = p + ggtitle("MDS on various distance metrics for ML16S dataset")
p
#this doesn't tell much other than 10 diff [15,18,25] diff 31

ML16S_unweighted_unifrac<-UniFrac(ML16S_final_biom_w_tree_w_metadata_w_refseq, weighted=FALSE)

#standardize
total = median(sample_sums(ML16S_final_biom_w_metadata_w_refseq))
standf = function(x, t=total) round(t*(x/sum(x)))
ML16S_standardized = transform_sample_counts(ML16S_final_biom_w_metadata_w_refseq, standf)

#UNIFRAC



library("phyloseq","ggplot2","doParallel","foreach")
library("ape")
library("ggdendro")
library("gridExtra")
set.seed(711)
#ggplot theme
# theme_set(theme_bw())
# fontsize = 18L
# theme_update(axis.title.x = element_text(size=fontsize))
# theme_update(axis.title.y = element_text(size=fontsize))
# theme_update(plot.title = element_text(size=fontsize+2))

theme_set(theme_classic())


ML16S_final_biom<-import_biom('otu_table_mc2_w_tax_noChloroplast_2012only_noSingletons_updated_taxonomy_NOmetadata2.biom')
ML16S_final_biom_w_tree<-import_biom('otu_table_mc2_w_tax_noChloroplast_2012only_noSingletons_updated_taxonomy_NOmetadata2.biom', treefilename = 'FINAL_biom_filtered_rep_set.tre')
ML16S_final_biom_w_metadata<-import_biom('otu_table_mc2_w_tax_noChloroplast_2012only_noSingletons_w_metadata.biom')
ML16S_final_biom_w_tree_w_metadata_w_refseq<-import_biom('otu_table_mc2_w_tax_noChloroplast_2012only_noSingletons_updated_taxonomy_w_metadata2.biom', treefilename='FINAL_biom_filtered_rep_set.tre', refseqfilename='FINAL_biom_filtered_rep_set_edited.fna', parseFunction=parse_taxonomy_greengenes)

ML16S_aplha<-plot_richness(ML16S_final_biom_w_metadata_w_refseq, x = "Depth", measures=c("Observed","Chao1","Shannon"), color="Depth")

ML16S_aplha+
  theme(axis.text.x=element_text(angle=90),
        legend.position="none",  
        panel.background = element_rect(fill=NA, color="black"))+
  ggtitle("Alpha Diversity")

ML16S_aplha+ggtitle("Alpha Diversity")+
  theme(strip.text.x = element_text(size = 8), plot.title = element_text(vjust = 2))+xlab("Depth (m)")
ML16S_alpha2<-ML16S_aplha+ggtitle("Alpha Diversity")+
  theme(strip.text.x = element_text(size = 8), plot.title = element_text(vjust = 2))+xlab("Depth (m)")

ML16S_aplha3<-plot_richness(ML16S_final_biom_w_metadata, x = "Depth", color="Depth")+
theme(axis.text.x=element_text(angle=90),
      legend.position="none",  
      panel.background = element_rect(fill=NA, color="black"))+
  ggtitle("Alpha Diversity")+xlab("Depth (m)")
ML16S_aplha3
#Unifrac, Weighted

ML16S_unweighted_unifrac<-UniFrac(ML16S_final_biom_w_tree_w_metadata_w_refseq, weighted=FALSE, fast=TRUE, parallel=FALSE)
ML16S_weighted_unifrac<-UniFrac(ML16S_final_biom_w_tree_w_metadata_w_refseq, weighted=TRUE, fast=TRUE, parallel=FALSE)

ML16S.hclust.unweighted     <- hclust(ML16S_unweighted_unifrac, method="average")
ML16S.hclust.weighted     <- hclust(ML16S_weighted_unifrac, method="average")
ML16S.hclust.unweighted$labels<-c("15 m", "18 m", "25 m", "10 m", "31 m")
ML16S.hclust.weighted$labels<-c("15 m", "18 m", "25 m", "10 m", "31 m")
ML16S.hclust.unweighted.dendr <- dendro_data(ML16S.hclust.unweighted, type="rectangle") 
ML16S.hclust.weighted.dendr<- dendro_data(ML16S.hclust.weighted, type="rectangle") 

ML16S.UF.UW<-ggplot() + 
  geom_segment(data=segment(ML16S.hclust.unweighted.dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(ML16S.hclust.unweighted.dendr), aes(x=x, y=y, label=label, hjust=0), size=3, hjust=-0.25) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())+
  ggtitle("Unweighted Unifrac Distance")
ML16S.UF.UW

ML16S.UF.W<-ggplot() + 
  geom_segment(data=segment(ML16S.hclust.weighted.dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(ML16S.hclust.weighted.dendr), aes(x=x, y=y, label=label, hjust=0), size=3, hjust=-0.25) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())+
  ggtitle("Weighted Unifrac Distance")
ML16S.UF.W

grid.arrange(ML16S_alpha2, arrangeGrob(ML16S.UF.UW, ML16S.UF.W,ncol=2), nrow=2)


setEPS()
postscript("Alpha_Beta.eps", paper="special",width=8, height=14, horizontal=FALSE)
showtext.begin()
grid.arrange(ML16S_alpha2, arrangeGrob(ML16S.UF.UW, ML16S.UF.W,ncol=2), nrow=2)
dev.off()


#below is old




labs = paste("sta_",1:50,sep="") #new labels
rownames(USArrests)<-labs #set new row names
hc <- hclust(dist(USArrests), "ave")



#convert cluster object to use with ggplot
dendr <- dendro_data(hc, type="rectangle") 

#your own labels (now rownames) are supplied in geom_text() and label=label




colorScale    <- rainbow(length(levels(get_variable(ML16S_final_biom_w_tree_w_metadata_w_refseq, "Depth"))))
cols          <- colorScale[get_variable(ML16S_final_biom_w_tree_w_metadata_w_refseq, "Depth")] 
tip.labels <- as(get_variable(ML16S_final_biom_w_tree_w_metadata_w_refseq, "Depth"), "character")

layout(matrix(c(1,2),1,2))
#plot WEIGHTED
plot(as.phylo(ML16S.hclust.weighted), show.tip.label = TRUE, tip.color = "white")
tiplabels(tip.labels, col = cols, frame = "none", adj = -0.05, cex=0.7)
#plot UNWEIGHTED
plot(as.phylo(ML16S.hclust.unweighted),show.tip.label=FALSE)
tiplabels(tip.labels, col = cols, frame = "none", adj = -0.05, cex=0.7)
xlab("Weighted UniFrac")



#tree (this is very messy - work on it?)
plot_tree(ML16S_final_biom_w_tree, size = "abundance")
#alpha diversity (current Fig. S1 as of 1-5-15)
plot_richness(ML16S_final_biom_w_metadata, x = "Depth")
#beta diversity
ML16S.weighted.unifrac<-UniFrac(ML16S_final_biom_w_tree, weighted=TRUE)
ML16S.unweighted.unifrac<-UniFrac(ML16S_final_biom_w_tree, weighted=FALSE)


