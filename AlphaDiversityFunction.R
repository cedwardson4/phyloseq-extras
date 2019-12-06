alphaDivSubsample<-function(ps,richness_measure="Observed",evenness_measure="InvSimpson",trials=100){
  #function implementing alpha diversity subsamped with replacement as described here:
  #http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#alpha_diversity
  #usage: ps = phyloseq object
  #richness_measure default is observed
  #evenness_measure default is inverse Simpson
  #available metrics from phyloseq see ?esimate_richness
  #n trials default=100
  library(phyloseq)
  min_lib<-min(sample_sums(ps))
  nsamp = nsamples(ps)
  trials = trials
  richness <- matrix(nrow = nsamp, ncol = trials)
  row.names(richness) <- sample_names(ps)
  evenness <- matrix(nrow = nsamp, ncol = trials)
  row.names(evenness) <- sample_names(ps)
  set.seed(100)

  for (i in 1:trials) {
    r <- rarefy_even_depth(ps, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
    # Calculate richness
   rich <- as.numeric(as.matrix(estimate_richness(r, measures = richness_measure)))
    richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = evenness_measure)))
  evenness[ ,i] <- even
  }

# Create a new dataframe to hold the means and standard deviations of richness estimates
  SampleID <- row.names(richness)
  mean <- apply(richness, 1, mean)
  sd <- apply(richness, 1, sd)
  measure <- rep("Richness", nsamp)
  rich_stats <- data.frame(SampleID, mean, sd, measure)

  # Create a new dataframe to hold the means and standard deviations of evenness estimates
  SampleID <- row.names(evenness)
  mean <- apply(evenness, 1, mean)
  sd <- apply(evenness, 1, sd)
  measure <- rep("Inverse Simpson", nsamp)
  even_stats <- data.frame(SampleID, mean, sd, measure)


  alpha <- rbind(rich_stats, even_stats)

  sample_data_df <- data.frame(sample_data(ps))
  alphadiv <- merge(alpha, sample_data_df, by = 1)
  return(alphadiv)
}

