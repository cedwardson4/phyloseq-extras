qiime2phyloseq <- function(otutablepath,treefilepath,mappingfilepath) {
  require(phyloseq)
  ps<-import_biom(otutablepath, 
                         treefilename=treefilepath,
                         parseFunction=parse_taxonomy_default)

  metadata<-import_qiime_sample_data(mappingfilepath)
  ps_combined<-merge_phyloseq(ps, metadata)

  colnames(tax_table(ps_combined))

  colnames(tax_table(ps_combined)) <- c(Rank1 = "Kingdom", Rank2 = "Phylum", Rank3 = "Class", 
                                           Rank4 = "Order", Rank5 = "Family", Rank6 = "Genus", Rank7 = "Species")

  tax_table(ps_combined)<-gsub("lca_D_0__", "lca_k_", tax_table(ps_combined))
  tax_table(ps_combined)<-gsub("lca_D_1__", "lca_p_", tax_table(ps_combined))
  tax_table(ps_combined)<-gsub("lca_D_2__", "lca_c_", tax_table(ps_combined))
  tax_table(ps_combined)<-gsub("lca_D_3__", "lca_o_", tax_table(ps_combined))
  tax_table(ps_combined)<-gsub("lca_D_4__", "lca_f_", tax_table(ps_combined))
  tax_table(ps_combined)<-gsub("lca_D_5__", "lca_g_", tax_table(ps_combined))
  tax_table(ps_combined)<-gsub("lca_D_6__", "lca_s_", tax_table(ps_combined))

  tax_table(ps_combined)<-gsub("D_[0-9]__", "", tax_table(ps_combined))
  return(ps_combined)
}
