alphaDivSubsample<-function(ps,richness_measure,evenness_measure,trials){
  min_lib<-min(sample_sums(ps))
  nsamp = nsamples(ps)
  trials = trials
  richness <- matrix(nrow = nsamp, ncol = trials)
  row.names(richness) <- sample_names(ps)
  evenness <- matrix(nrow = nsamp, ncol = trials)
  row.names(evenness) <- sample_names(ps)
  set.seed(100)

  for (i in 1:100) {
    r <- rarefy_even_depth(ps, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
    # Calculate richness
   rich <- as.numeric(as.matrix(estimate_richness(r, measures = richness_measure)))
    richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = evenness_measure)))
  evenness[ ,i] <- even
  }

# Create a new dataframe to hold the means and standard deviations of richness estimates
  X.SampleID <- row.names(richness)
  mean <- apply(richness, 1, mean)
  sd <- apply(richness, 1, sd)
  measure <- rep("Richness", nsamp)
  rich_stats <- data.frame(X.SampleID, mean, sd, measure)

  # Create a new dataframe to hold the means and standard deviations of evenness estimates
  X.SampleID <- row.names(evenness)
  mean <- apply(evenness, 1, mean)
  sd <- apply(evenness, 1, sd)
  measure <- rep("Inverse Simpson", nsamp)
  even_stats <- data.frame(X.SampleID, mean, sd, measure)


  alpha <- rbind(rich_stats, even_stats)

  sample_data_df <- data.frame(sample_data(ps))
  #AQComp_sample_df<-plyr::rename(AQComp_sample_df, replace=c("X.SampleID"="X.SampleID"))
  alphadiv <- merge(alpha, sample_data_df, by = "X.SampleID")
  return(alphadiv)
}

