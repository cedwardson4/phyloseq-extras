create_phyloseq_other_taxa<-function(physeq){
  #take a phyloseq object, separate into individual parts, add "Other" OTU row in tax_table, otu_table) 
  #physeq_with_other<-create_phyloseq_other_taxa(physeq)
  library(phyloseq)
  otu_table<-as(phyloseq::otu_table(physeq,taxa_are_rows=phyloseq::taxa_are_rows(phyloseq::otu_table(physeq))),"matrix")
  tax_table<-as(phyloseq::tax_table(physeq),"matrix")

  taxa_are_rows_logical<-phyloseq::taxa_are_rows(phyloseq::otu_table(physeq))
  
  #add "other" row, as difference in sum of rows, to otu_table
  if(taxa_are_rows_logical==TRUE) {
    otu_table_w_other<-rbind(otu_table,apply(otu_table, 2,function(x) 1-sum(x)))
    rownames(otu_table_w_other)[nrow(otu_table_w_other)]<-"ZOther"
  } else if (taxa_are_rows_logical==FALSE) {
    otu_table_w_other<-cbind(otu_table,apply(otu_table, 1,function(x) 1-sum(x)))
    colnames(otu_table_w_other)[ncol(otu_table_w_other)]<-"ZOther"
  }
  #add "other" row to tax_table
  other_tax_mat<-matrix(replicate(ncol(tax_table),"ZOther"),nrow=1,ncol=ncol(tax_table))
  tax_table_w_other<-rbind(tax_table,other_tax_mat)
  rownames(tax_table_w_other)[nrow(tax_table_w_other)]<-"ZOther"
  #merge and return new object
  physeq_with_other<-phyloseq::merge_phyloseq(phyloseq::otu_table(otu_table_w_other,taxa_are_rows=taxa_are_rows_logical),
                                              phyloseq::tax_table(tax_table_w_other),sample_data(physeq))
  return(physeq_with_other)
}