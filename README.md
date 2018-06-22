# phyloseq-extras
Scripts that I wrote or borrowed from others that I find useful for doing stuff with phyloseq objects(or more generally) for microbiome analyses.

## phyloseq2qiime2
Requries the following R packages (and loads them):

[phyloseq](https://github.com/joey711/phyloseq)

[ape](https://cran.r-project.org/web/packages/ape/index.html)

[biomformat](https://bioconductor.org/packages/release/bioc/html/biomformat.html) > 1.7

[Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)

This script takes a phyloseq object, determines if each of the slots (otu_table, tax_table, sam_data, phy_tree, refseqs) are present, and if so, writes them to a format that is importable by qiime2.

## qiime2phyloseq
This script takes a biom, tree, and, mapping file processed by my QIIME 1 pipeline and automatically imports and formats a phyloseq object. It probably won't work for anyone else right now. Adding it here for when I get my QIIME 1 pipeline up here or to generalize. 

## Phyloseq_Top_Percent (function is phyloseq_top_percent_taxa())
This function takes a normalized (normalize=FALSE, default), or unnormalized (normalize=TRUE) OTU table and subsets only taxa above some percent (0-100%) of taxa that is found in n samples (kOverA function).
Example: 
If you want all taxa present in one or more samples >2%
```
ps_min2percent<-phyloseq_top_percent_taxa(ps, n_min=1, percent=2))
```
This will also print out the number of taxa at the OTU/ASV, Genus and Family level to indicate number of colors needed for scale_fill_manual with custom palettes.

## Phyloseq_Create_Other (function is create_phyloseq_other_taxa())
This function takes a subsampled normaizled (relative abundance) otu_table and tax_table from a phyloseq object and adds an "other" taxa category so all samples sum to 1.
- Issues: currently this adds the "Other" taxa as "ZOther" so that it appears at the bottom of the legend (unless you have Zymomonas or other Z taxa!). I need to add factor leveling to make this cleaner and work properly. 

## plot_ordered_bar2
The second example of ordered plot_bar from [here](https://github.com/joey711/phyloseq/issues/442)

## AlphaDiversityFunction (function is alphaDivSubsample())
Takes a phyloseq object and performs richness and eveness calculations using subsampling as described [here](http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#alpha_diversity).

## veganotu
From https://github.com/joey711/phyloseq/issues/190. Make phyloseq otu_table usable in general vegan functions.

## core_microbiome
From https://github.com/davidelliott/core-microbiome/

## function_bio.env & function_bv.step
BIOENV and BV-STEP for "comparison of distance/similarity matrices between two sets of data having either samples or variables in common" 
From http://menugget.blogspot.com/2011/06/clarke-and-ainsworths-bioenv-and-bvstep.html

## import_biom2
Import biom files to phyloseq.
See https://gist.github.com/jnpaulson/324ac1fa3eab1bc7f845

## miseqR
Potentially useful functions from the Denef lab
https://github.com/michberr/MicrobeMiseq

## read_hdf5_biom
Imports an hdf5-formatted BIOM file. From https://gist.github.com/jnpaulson/c9cdcb8cf95e7daae468. This may not be necessary as it appears the r package biomformat now has a read_hdf5_biom function.

## runDada2
Likely outdated and incorrect automated dada2 pipeline. Might be useful if updated for current dada2.
From https://github.com/zachcp/phyloseq-tools.

## VarDistance (function varDistance.f())
From [Ashley Shade](https://github.com/ashleyshade). Not specifcally for phyloseq objects, but useful for beta diversity.

## venn_diagram2
From https://github.com/hallamlab/mp_tutorial/blob/master/downstream_analysis_r/code/venn_diagram2.r
Another non-phyloseq specific function, but may be useful for making venn diagrams.
