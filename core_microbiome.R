#https://github.com/davidelliott/core-microbiome/blob/master/core-microbiome.Rmd

# core microbiome function
# finds the core microbiome with various options
# physeq: the phyloseq object to query
# threshold: the abundance below which OTUs are not regarded as being "core taxa"
# method: "max" or "mean" - define whether to calculate the threshold for the max or mean
# ranks: vector of taxonomic ranks at which to search for the core taxa
# group: factor on which to merge samples before applying threshold criteria. 

# Note.
# max and mean are calculated for the whole OTU table per row
# therefore
# mean gives the overall abundant OTUs
# max will include any OTU which is abundant any 1 sample, so it is really like a sample-wise search

# returns a phyloseq object containing only the OTUs matching the criteria set
core_microbiome <- function(physeq,threshold=2,method="mean",ranks=c("Species"),group="",verbose=0){
  core_otus_glom <- c()
  physeq.original <- physeq
  if(group!="") {
    # note that grouping will sum the counts
    # but for the threshold to be applied correctly we need to
    # use the mean for the group.
    # the mean is worked out and applied below.
    # I'm sure this can be done in a better way
    # HAS NOT BEEN TESTED ON OBJECTS WHERE THE GROUP FREQUENCY VARIES
    sample_data(physeq)$group_frequency <- NA
    for(n in 1:nrow(sample_data(physeq))) {
      sample_data(physeq)$group_frequency <- sum(get_variable(physeq, group)==get_variable(physeq, group)[n])
    }
    otu_table(physeq) <- otu_table(physeq)/sample_data(physeq)$group_frequency
    physeq <- merge_samples(x=physeq,group=group)
    # I don't know why the above transposes the OTU table. Turn it back:
    otu_table(physeq) <- t(otu_table(physeq))
    
    # transform_sample_counts(physeq, function(x) x/sample_data(physeq)$group_frequency)
  }
  for (rank in ranks) {
    # agglomerate to the taxonomic rank being searched
    # but only do this if we are not at the bottom rank already
    if(rank!=rank_names(physeq)[length(rank_names(physeq))]) {
      physeq.glom <- tax_glom(physeq,rank)
    } else {
      physeq.glom <- physeq
    }
    # get the OTU table into a dataframe
    glom_otus <- as.data.frame(otu_table(physeq.glom))
    # calculate the test statistic for each OTU
    glom_otus$stat <- apply(glom_otus,1,method)
    # make a list of the taxa meeting the abundance criteria
    glom_rank_candidates <- row.names(glom_otus[glom_otus$stat >= threshold,])
    # identify the taxonomy of the candidates
    if(length(glom_rank_candidates)>0) {
      glom_rank_taxonomy <- tax_table(physeq.glom)[glom_rank_candidates,rank]
      
      # identify the most abundant OTU within the identified taxonomic group, and add it to the list
      # tax_glom does not always pick the most abundant OTU in the classification
      # Therefore need to specifically check for the most abundant OTU 
      # should only do this if we are not working at the lowest rank in the phyloseq object
      if(rank!=rank_names(physeq)[length(rank_names(physeq))]) {
        for (candidate in glom_rank_taxonomy) {
          OTUs_in_taxonomic_group <- tax_table(physeq)[,rank] == as.character(candidate)
          most_abundant_OTU <- names(which.max(taxa_sums(physeq)[as.vector(OTUs_in_taxonomic_group)]))
          core_otus_glom <- c(core_otus_glom,most_abundant_OTU)
        }
      } 
      else {
        # at OTU level we don't need to search for the most abundant one - we want them all
        core_otus_glom <- c(core_otus_glom,glom_rank_candidates)  
      }
    }
    if(verbose>1)  print(paste(rank,"level search found",length(glom_rank_candidates),"taxa above",method,"abundance of",threshold))
  }
  if(verbose>0) print(paste("Search found", length(core_otus_glom),"unique OTUs"))
  core_otus <- unique(core_otus_glom)
  return(prune_taxa(core_otus,physeq.original) )
}


# core compare function
# Input 2 phyloseq objects which have been pruned to contain a 'core microbiome'
# returns a list containing:
# extra OTUs in the first object, not present in the second
# extra taxonomic groups represented by those OTUs, which are not present in the second object

core_compare <- function(physeq1,physeq2) {
  out <- list()
  # Find out which OTUs were added by the multi-rank approach
  extras <- tax_table(physeq1)[!rownames(otu_table(physeq1))%in%rownames(otu_table(physeq2))]
  out[["extra_otus"]] <- rownames(extras)
  # Find out which taxa are additionally represented by those OTUs
  for(rank in rank_names(physeq1)) {
    index <- !extras[,rank]%in%unique(tax_table(physeq2)[,rank])
    if(sum(index)) {
      out[[rank]] <- unique(as.vector(extras[index,rank]))
    }
    else {
      out[rank] <- "none"
    }
  }
  return(out)
}
