# phyloseq-extras
Scripts for doing stuff with phyloseq objects.

## phyloseq2qiime2
Requries the following R packages (and loads them):

[phyloseq](https://github.com/joey711/phyloseq)

[ape](https://cran.r-project.org/web/packages/ape/index.html)

[biomformat](https://bioconductor.org/packages/release/bioc/html/biomformat.html) > 1.7

[Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)

This script takes a phyloseq object, determines if each of the slots (otu_table, tax_table, sam_data, phy_tree, refseqs) are present, and if so, writes them to a format that is importable by qiime2.

## AlphaDiversityFunction (function is alphaDivSubsample())
Takes a phyloseq object and performs richness and eveness calculations using subsampling as described [here](http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#alpha_diversity).
