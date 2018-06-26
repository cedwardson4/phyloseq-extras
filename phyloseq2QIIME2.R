phyloseq2qiime2<-function(physeq){
  #take a phyloseq object,check for individual parts, write to files ready for qiime2 upload
  library(phyloseq)
  library(biomformat)
  library(ape)
  library(Biostrings)
  if(packageVersion("biomformat") < "1.7") {
    stop("This will only work with biomformat version > 1.7")
  }
  ps_name <-deparse(substitute(physeq))
  taxa_are_rows_logical<-taxa_are_rows(physeq)
  #write OTU table to biom file
  if(is.null(access(physeq,"otu_table"))==FALSE){
    if(taxa_are_rows_logical==TRUE) {
      otu<-as(otu_table(physeq),"matrix")
      otu_biom<-make_biom(data=otu)
      write_biom(otu_biom,biom_file=paste0(ps_name,"_features-table.biom"))
      print(paste0("Writing feature table to ",ps_name,"_feature-table.biom"))
    } else if (taxa_are_rows_logical==FALSE) {
      otu<-t(as(otu_table(physeq),"matrix"))
      otu_biom<-make_biom(data=otu)
      write_biom(otu_biom,biom_file=paste0(ps_name,"_feature-table.biom"))
      print(paste0("Writing feature table to ",ps_name,"_feature-table.biom"))
    }
  }
  #write sample data (metadata) to tsv
  if(is.null(access(physeq,"sam_data"))==FALSE){
    write.table(sample_data(physeq),file=paste0(ps_name,"_sample-metadata.txt"), 
                sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    print(paste0("Writing sample metadata to ",ps_name,"_sample-metadata.txt"))
  }
  #write taxonomy table to qiime2 formatted taxonomy
  if(is.null(access(physeq,"tax_table"))==FALSE){
    tax<-as(tax_table(physeq),"matrix")
    tax_cols <- colnames(tax)
    tax<-as.data.frame(tax)
    tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
    for(co in tax_cols) tax[co]<-NULL
    write.table(tax, file=paste0(ps_name,"_tax.txt"), quote=FALSE, col.names=FALSE, sep="\t")
    print(paste0("Writing taxonomy table to ",ps_name,"_tax.txt"))
  }
  #write phylogenetic tree to newick formwat
  if(is.null(access(physeq,"phy_tree"))==FALSE){
    if(is.rooted(phy_tree(physeq))==TRUE) {
      ape::write.tree(phy_tree(physeq),file=paste0(ps_name,"_tree-rooted.newick"))
      print(paste0("Writing rooted tree to ",ps_name,"_tree-rooted.newick"))
    } else if (is.rooted(phy_tree(physeq))==FALSE) {
      ape::write.tree(phy_tree(physeq),file=paste0(ps_name,"_tree-unrooted.newick"))
      print(paste0("Writing unrooted tree to ",ps_name,"_tree-unrooted.newick"))
    }      
  }
  #write representative sequences to fasta format
  if(is.null(access(physeq,"refseq"))==FALSE){
    writeXStringSet(refseq(physeq),filepath=paste0(ps_name,"_ref-seqs.fasta"))
    print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))
  } else if (taxa_are_rows_logical==FALSE && unique(grepl("[^ATCG]",colnames(otu_table(physeq)),ignore.case=TRUE) == FALSE)) {
    uniquesToFasta(t(otu), fout=paste0(ps_name,"_ref-seqs.fasta"), ids=rownames(otu))
    print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))      
  } else if (taxa_are_rows_logical==TRUE && unique(grepl("[^ATCG]",rownames(otu_table(physeq)),ignore.case=TRUE) == FALSE)) {
    uniquesToFasta(otu, fout=paste0(ps_name,"_ref-seqs.fasta"), ids=rownames(otu))
    print(paste0("Writing reference sequences to FASTA file ",ps_name,"_ref-seqs.fasta"))    
  }
}
