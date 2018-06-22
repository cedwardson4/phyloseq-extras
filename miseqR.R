##### FUNCTIONS ########
#This is a set of functions for analyzing microbial amplicon data written by members of the Denef lab
#https://github.com/michberr/MicrobeMiseq



##### Normalization #######

# Scales reads by 
# 1) taking proportions,
# 2) multiplying by a given library size of n
# 3) rounding down
scale_reads <- function(physeq, n) {
  physeq.scale <-
    transform_sample_counts(physeq, function(x) {
      (n * x/sum(x))
    })
  otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}

##### ADONIS ###########

# Function to run adonis test on a physeq object and a variable from metadata 
doadonis <- function(physeq, category) {
  physeq.scale <- scale_reads(physeq, min(sample_sums(physeq)))
  bdist <- phyloseq::distance(physeq.scale, "bray")
  col <- as(sample_data(physeq), "data.frame")[ ,category]
  
  # Adonis test
  adonis.bdist <- adonis(bdist ~ col)
  print("Adonis results:")
  print(adonis.bdist)
  
  # Homogeneity of dispersion test
  betatax = betadisper(bdist,col)
  p = permutest(betatax)
  print("Betadisper results:")
  print(p$tab)
}


########## Bar Plots #################

# This function takes a phyloseq object, agglomerates OTUs to the desired taxonomic rank, 
# prunes out OTUs below a certain relative proportion in a sample (ie 1% ) 
# and melts the phyloseq object into long format.
transform_and_melt <- function(physeq, taxrank, prune) {
  
  # Agglomerate all otu's by given taxonomic level
  physeq_taxrank <- tax_glom(physeq, taxrank = taxrank)
  
  # Create a new phyloseq object which removes taxa from each sample below the prune parameter
  physeq_taxrank.prune <- transform_sample_counts(physeq_taxrank,function(x) {x/sum(x)})
  otu_table(physeq_taxrank.prune)[otu_table(physeq_taxrank.prune) < prune] <- 0
  physeq_taxrank.prune <- prune_taxa(taxa_sums(physeq_taxrank.prune) > 0, physeq_taxrank.prune)
  
  # Melt into long format and sort by sample and taxrank
  physeq_taxrank.long <- psmelt(physeq_taxrank.prune)
  names(physeq_taxrank.long)[names(physeq_taxrank.long) == taxrank] <- "Taxonomy"
  physeq_taxrank.long <- arrange(physeq_taxrank.long, Sample, Taxonomy)
  
  # Return long data frame
  return(physeq_taxrank.long)
}


# This function takes  a data frame in long format 
# (such as the output from transform_and_melt) 
# and produces a stacked barplot of community composition.
make_tax_barplot <- function(df, x, y, tax, facet, title, colors, xlab, ylab, relative, outline, guide) {
  ggplot(df, aes_string(x = x, y = y, fill = tax)) +
    facet_grid(reformulate(facet), scales="free_y") +
    geom_bar(stat = "identity", show_guide = guide) +
    scale_fill_manual(values = colors) +
    scale_x_discrete(
      breaks = c("6/10", "7/8", "8/4", "9/2", "10/6", "11/3"),
      labels = c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov"), drop = FALSE
    ) +
    theme(
      axis.title.x = element_text(size = 16,face = "bold"),
      axis.text.x = element_text(angle = 50, colour = "black", vjust = 1, hjust = 1, size = 13),
      axis.text.y = element_text(colour = "black", size = 10),
      axis.title.y = element_text(face = "bold", size = 16),
      plot.title = element_text(size = 18),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      legend.position = "right",
      strip.text.x = element_text(size = 16, face = "bold"),
      strip.text.y = element_text(size = 16, face = "bold"),
      strip.background = element_rect(color = "white",size = 2, fill = NA),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
      panel.margin = unit(1, "lines"),
      panel.background = element_rect(fill = "#a8a8a8")
    ) +
    guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(title) +
    if (relative) {
      geom_bar(
        stat = "identity",position = "fill", colour = outline,show_guide = FALSE
      )
    } else {
      geom_bar(stat = "identity",colour = outline,show_guide = FALSE)
    }
  
}



############## Ordinations ######################

# Wrapper function for creating an ordination
# 1) scales an OTU table 
# 2) Ordinates
# 3) Plots ordination

myord <- function(physeq, n,
                  method, distance, colors = NULL, 
                  factor.color, factor.shape, title){
  
  physeq.scale <- scale_reads(physeq, n)
  physeq.ordinate <- ordinate(physeq.scale, method = method, distance = distance)
  
  # Plot the ordination
  myordplot(physeq.scale, physeq.ordinate, colors, factor.color, factor.shape, title)
}


myordplot <- function (physeq, ordination, colors, factor.color, factor.shape, title){
  plot_ordination(
    physeq = physeq, 
    ordination = ordination, 
    color = factor.color, 
    shape = factor.shape) +
    geom_point(
      aes_string(color = factor.color), 
      alpha = 0.7, 
      size = 4) +
    geom_point(
      colour="grey90", 
      size = 1.5) + 
    ggtitle(title) +
    if(!is.null(colors)){
      scale_color_manual(values = colors)
    }
  
}

###### Merge functions ############

#Merge samples by averaging OTU countsinstead of summing
merge_samples_mean <- function(physeq, group){
  # Calculate the number of samples in each group
  group_sums <- as.matrix(table(sample_data(physeq)[ ,group]))[,1]
  
  # Merge samples by summing
  merged <- merge_samples(physeq, group)
  
  # Divide summed OTU counts by number of samples in each group to get mean
  # Calculation is done while taxa are columns
  x <- as.matrix(otu_table(merged))
  if(taxa_are_rows(merged)){ x<-t(x) }
  out <- floor(t(x/group_sums))
  
  # Return new phyloseq object with
  
  out <- otu_table(out, taxa_are_rows = TRUE)
  otu_table(merged) <- out
  return(merged)
}

# Merge samples, just including OTUs that were present in all

masa_merge <- function(physeq, group){
  
  # Make sure we're not starting with more taxa than we need 
  physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
  
  s <- data.frame(sample_data(physeq))
  l <- levels(s[,group])
  o <- otu_table(physeq)
  
  # Loop through each category
  for (cat in 1:length(l)) {
    
    # Get the index of all samples in that category
    w <- which(s[,group]==l[cat])
    
    # subset to just those columns of OTU table
    cat.sub<-o[,w]
    print(dim(cat.sub))
    
    # Find the indices of 0's in the OTU table
    zeros <- apply(cat.sub, 1, function(r) any(r==0))
    
    # If an OTU had a 0 in at least one sample, change all samples to 0
    cat.sub[zeros,] <- 0
  }
  
  o[,w] <- cat.sub
  otu_table(physeq) <- o
  
  return(physeq)
  
}



