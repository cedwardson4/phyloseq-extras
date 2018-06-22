import_biom2 <- function(x, 
                         treefilename=NULL, refseqfilename=NULL, refseqFunction=readDNAStringSet, refseqArgs=NULL,
                         parseFunction=parse_taxonomy_default, parallel=FALSE, version=1.0, ...){
  
  # initialize the argument-list for phyloseq. Start empty.
  argumentlist <- list()
  
  x = biom(x)
  b_data = biom_data(x)
  b_data_mat = as(b_data, "matrix")
  otutab = otu_table(b_data_mat, taxa_are_rows=TRUE)
  argumentlist <- c(argumentlist, list(otutab))
  
  ########################################
  # Taxonomy Table
  ########################################
  # Need to check if taxonomy information is empty (minimal BIOM file)
  if(  all( sapply(sapply(x$rows, function(i){i$metadata}), is.null) )  ){
    taxtab <- NULL
  } else {
    # parse once each character vector, save as a list
    taxlist = lapply(x$rows, function(i){
      parseFunction(i$metadata$taxonomy)
    })
    names(taxlist) = sapply(x$rows, function(i){i$id})
    taxtab = build_tax_table(taxlist)
  }
  argumentlist <- c(argumentlist, list(taxtab))
  
  ########################################
  # Sample Data ("columns" in QIIME/BIOM)
  ########################################
  # If there is no metadata (all NULL), then set sam_data <- NULL
  if( is.null(sample_metadata(x)) ){
    samdata <- NULL
  } else {
    samdata = sample_data(sample_metadata(x))
  }
  argumentlist <- c(argumentlist, list(samdata))
  
  ########################################
  # Tree data
  ########################################
  if( !is.null(treefilename) ){
    if( inherits(treefilename, "phylo") ){
      # If argument is already a tree, don't read, just assign.
      tree = treefilename
    } else {
      # NULL is silently returned if tree is not read properly.
      tree <- read_tree(treefilename, ...)		
    }
    # Add to argument list or warn
    if( is.null(tree) ){
      warning("treefilename failed import. It not included.")
    } else {
      argumentlist <- c(argumentlist, list(tree) )
    }
  }
  
  ########################################
  # Reference Sequence data
  ########################################	
  if( !is.null(refseqfilename) ){
    if( inherits(refseqfilename, "XStringSet") ){
      # If argument is already a XStringSet, don't read, just assign.
      refseq = refseqfilename
    } else {
      # call refseqFunction and read refseqfilename, either with or without additional args
      if( !is.null(refseqArgs) ){
        refseq = do.call("refseqFunction", c(list(refseqfilename), refseqArgs))
      } else {
        refseq = refseqFunction(refseqfilename)
      }
    }
    argumentlist <- c(argumentlist, list(refseq) )
  }
  
  ########################################
  # Put together into a phyloseq object
  ########################################
  return( do.call("phyloseq", argumentlist) )
  
}