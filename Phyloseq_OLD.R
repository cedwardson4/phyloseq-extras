setwd("F:/R")
#if(!exists("import_biom2", mode="function")) source("import_biom2.R")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library("dplyr")
library("ggdendro")
library("scales")
library("grid")
library("vegan")
library("reshape2")
library("dtplyr")


#library("biom"); packageVersion("biom")
#library("rhdf5"); packageVersion("rhdf5")
theme_set(theme_bw())

YR250PE_OTUs<-import_biom("YR250-PE/otu_table_mc2_w_tax_no_pynast_failures_no_cmu.biom", 
                          treefilename='YR250-PE/rep_set_mc2_no_pynast_failures_no_cmu.tre', 
                          parseFunction=parse_taxonomy_default)

YR250PE_OTUs_map<-import_qiime_sample_data("YR250-PE/YR-Run250PE-MappingFile.txt")

YR250PE_OTUs_meta<-merge_phyloseq(YR250PE_OTUs, YR250PE_OTUs_map)
colnames(tax_table(YR250PE_OTUs_meta)) 

colnames(tax_table(YR250PE_OTUs_meta)) <- c(Rank1 = "Kingdom", Rank2 = "Phylum", Rank3 = "Class", 
                                  Rank4 = "Order", Rank5 = "Family", Rank6 = "Genus", Rank7 = "Species")

tax_table(YR250PE_OTUs_meta)<-gsub("D_[0-9]__", "", tax_table(YR250PE_OTUs_meta))
head(tax_table(YR250PE_OTUs_meta))

YR250PE_OTUs_meta

#Phyloseq Object From DADA2 pipeline
YR250PE_ps


#read count distribution
sample_sum_df_qiime <- data.frame(sum = sample_sums(YR250PE_OTUs_meta))
sample_sum_df_dada2 <- data.frame(sum = sample_sums(YR250PE_ps))

# Histogram of sample read counts 1
ggplot(sample_sum_df_qiime, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

# Histogram of sample read counts 2

sdt = data.table(as(sample_data(YR250PE_ps), "data.frame"),
                 TotalReads = sample_sums(YR250PE_ps), keep.rownames = TRUE)

YR_ps_seq_depth<-ggplot(sdt, aes(TotalReads)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
YR_ps_seq_depth
YR_ps_seq_depth+facet_wrap(~BodySite)


ggplot(sdt, aes(SamplingDate, TotalReads, color = as.factor(SamplingLocation))) + 
  geom_point(size = 5) + 
  geom_smooth(method = lm) +
  ggtitle("Sequencing Depth vs. Time") +
  scale_y_log10() +
  facet_grid(BodySite ~ .)


# mean, max and min of sample read counts
smin <- min(sample_sums(YR250PE_ps))
smean <- mean(sample_sums(YR250PE_ps))
smax <- max(sample_sums(YR250PE_ps))

# melt to long format (for ggploting) 
# prune out phyla below 2% in each sample

YR250PE_ps_phylum <- YR250PE_ps %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

YR250PE_OTUs_meta_phylum <- YR250PE_OTUs_meta %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)


ggplot(YR250PE_ps_phylum, aes(x = BodySite, y = Abundance, fill = Phylum)) + 
  facet_grid(Caught_Born_Habitat~.) +
  geom_bar(stat = "identity", position="fill") +
  scale_fill_manual(values = phylum_colors) +
  # scale_x_discrete(
  #   breaks = c("7/8", "8/4", "9/2", "10/6"),
  #   labels = c("Jul", "Aug", "Sep", "Oct"), 
  #   drop = FALSE
  # ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition of Yellow Rays \n Bacterial Communities by Sampling Site") 


#Unconstrained Ordinations
# Scale reads to even depth 
#YR250PE_ps_scale <- scale_reads(YR250PE_ps)
#code didn't work...use builtin

YR250PE_ps_even10k <- rarefy_even_depth(YR250PE_ps, rngseed=TRUE, sample.size=10000, replace=FALSE)

# Ordinate
YR250PE_ps_pcoa <- ordinate(
  physeq = YR250PE_ps_even10k, 
  method = "PCoA", 
  distance = "bray"
)

# Plot 
plot_ordination(
  physeq = YR250PE_ps_even10k,
  ordination = YR250PE_ps_pcoa,
  color = "SamplingLocation",
  shape = "BodySite",
  title = "PCoA of Yellow Ray bacterial Communities"
) + 
  #scale_color_manual(values = c("#a65628", "red", "#ffae19",
#                                "#4daf4a", "#1919ff", "darkorchid3", "magenta")
  #) +
  geom_point(aes(color = SamplingLocation), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 


#NMDS
set.seed(1)
YR250PE_ps_nmds <- ordinate(
  physeq = YR250PE_ps_even10k, 
  method = "NMDS", 
  distance = "bray"
)

plot_ordination(
  physeq = YR250PE_ps_even10k,
  ordination = YR250PE_ps_nmds,
  color = "SamplingLocation",
  shape = "BodySite",
  title = "NMDS of Yellow Ray bacterial Communities"
) + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19",
                                "#4daf4a", "#1919ff", "darkorchid3", "magenta")
  ) +
  geom_point(aes(color = SamplingLocation), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 


#permanova

set.seed(1)

# Calculate bray curtis distance matrix
YR250PE_ps_bray <- phyloseq::distance(YR250PE_ps_even10k, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(YR250PE_ps_even10k))

# Adonis test
adonis(YR250PE_ps_bray ~ SamplingLocation+BodySite+Caught_Born+Habitat, data = sampledf)

# Homogeneity of dispersion test
beta <- betadisper(YR250PE_ps_bray, sampledf$SamplingLocation)
permutest(beta)

#Constrained ordination
# CAP ordinate
cap_ord <- ordinate(
  physeq = YR250PE_ps_even10k, 
  method = "CAP",
  distance = YR250PE_ps_bray,
  formula = ~ SamplingLocation + BodySite + Habitat + Caught_Born_Habitat
)

# CAP plot
cap_plot <- plot_ordination(
  physeq = YR250PE_ps_even10k, 
  ordination = cap_ord, 
  color = "SamplingLocation", 
  axes = c(1,2)
) + 
  aes(shape = BodySite) + 
  geom_point(aes(colour = SamplingLocation), alpha = 0.4, size = 4) + 
  geom_point(colour = "grey90", size = 1.5) + 
  scale_color_manual(values = c("#a65628", "red", "#ffae19", "#4daf4a", 
                                "#1919ff", "darkorchid3", "magenta")
  )


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )

anova(cap_ord)

###Phyloseq ordination






#calculate distance metrics (run all or subset)
#remove DPCoA and ANY
dist_methods <- dist_methods[-(3:3)]
dist_methods <- dist_methods[-which(dist_methods=="ANY")]
dist_methods <- dist_methods[-which(dist_methods=="manhattan")]

YR250PE_distance <- vector("list", length(dist_methods))
names(YR250PE_distance) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- distance(YR250PE_OTUs_meta, method=i)
  # Calculate ordination
  iMDS  <- ordinate(YR250PE_OTUs_meta, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(YR250PE_OTUs_meta, iMDS, color="Caught_Born_Habitat", shape="BodySite")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  YR250PE_distance[[i]] = p
}

df = ldply(YR250PE_distance, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=Caught_Born_Habitat, shape=BodySite))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for Yellow Ray 250PE dataset")
p


#testDist <- distance(YR250PE_OTUs_meta, method="manhattan")
#testMDS  <- ordinate(YR250PE_OTUs_meta, "MDS", distance=testDist)


YR250PE_OTUs_UniFrac_unweighted<-UniFrac(YR250PE_OTUs_meta, weighted=FALSE)
YR250PE_OTUs_UniFrac_weighted<-UniFrac(YR250PE_OTUs_meta, weighted=TRUE)

YR250PE.hclust.unweighted     <- hclust(YR250PE_OTUs_UniFrac_unweighted, method="average")
YR250PE.hclust.weighted     <- hclust(YR250PE_OTUs_UniFrac_weighted, method="average")
YR250PE.hclust.unweighted.dendr <- dendro_data(YR250PE.hclust.unweighted, type="rectangle") 
YR250PE.hclust.weighted.dendr<- dendro_data(YR250PE.hclust.weighted, type="rectangle") 

YR250PE.UF.UW<-ggplot() + 
  geom_segment(data=segment(YR250PE.hclust.unweighted.dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(YR250PE.hclust.unweighted.dendr), aes(x=x, y=y, label=label, hjust=0), size=3, hjust=-0.25) +
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
YR250PE.UF.UW

#DESEQ2 through Phyloseq

library("DESeq2"); packageVersion("DESeq2")
YR250PE_OTUs_meta_subset <- subset_samples(YR250PE_OTUs_meta, Habitat != "None")
YR250PE_OTUs_meta_subset <- prune_samples(sample_sums(YR250PE_OTUs_meta) > 500, YR250PE_OTUs_meta)
YR250PE_OTUs_meta_subset
YR250PE_OTUs_deseq <- phyloseq_to_deseq2(YR250PE_OTUs_meta_subset, ~ Habitat)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}



geoMeans = apply(counts(YR250PE_OTUs_deseq), 1, gm_mean)
YR250PE_OTUs_deseq = estimateSizeFactors(YR250PE_OTUs_deseq, geoMeans = geoMeans, type="iterate")
YR250PE_OTUs_deseq = DESeq(YR250PE_OTUs_deseq, test="Wald", fitType="parametric")

rs <- rowSums(counts(YR250PE_OTUs_deseq))
rmx <- apply(counts(YR250PE_OTUs_deseq), 1, max)
plot(rs+1, rmx/rs, log="x")



res = results(YR250PE_OTUs_deseq, cooksCutoff = FALSE)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(YR250PE_OTUs_meta_subset)[rownames(sigtab), ], "matrix"))
head(sigtab)

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

negsigtab = sigtab[sigtab[, "log2FoldChange"] < 0, ]
negsigtab = negsigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Class))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
