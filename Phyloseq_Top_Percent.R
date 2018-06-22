phyloseq_top_percent_taxa<-function(physeq,n_min,percent,normalize=FALSE){
  #function takes percent, kOverA takes proportion
  library(phyloseq)
  library(genefilter)
  proportion=percent/100
  filter_function<-genefilter::filterfun(genefilter::kOverA(n_min,proportion))
  #if physeq is NOT relative abundance
  if(mean(phyloseq::sample_sums(physeq))!=1 & normalize==TRUE){
    physeq<-phyloseq::transform_sample_counts(physeq,function(x) x / sum (x))
    physeq_filter_list<-filter_taxa(physeq,filter_function)
    physeq_filtered<-prune_taxa(physeq_filter_list,physeq)
  } else if(mean(phyloseq::sample_sums(physeq))!=1 & normalize==FALSE) {
    physeq_filter_list<-filter_taxa(physeq,filter_function)
    physeq_filtered<-prune_taxa(physeq_filter_list,physeq)
  } else if( mean(phyloseq::sample_sums(physeq))==1) {
    physeq_filter_list<-filter_taxa(physeq,filter_function)
    physeq_filtered<-prune_taxa(physeq_filter_list,physeq)
  }
    print(paste("The number of colors needed for OTU coloring is", length(rownames(tax_table(physeq_filtered)))))
    print(paste("The number of colors needed for Genus coloring is", length(unique((tax_table(physeq_filtered)[,"Genus"])))))
    print(paste("The number of colors needed for Family coloring is", length(unique((tax_table(physeq_filtered)[,"Family"])))))
    
  return(physeq_filtered)
}

  