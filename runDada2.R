#!/usr/bin/env Rscript

#from: https://raw.githubusercontent.com/zachcp/phyloseq-tools/master/R/scripts/runDada2.R

suppressPackageStartupMessages(library(docopt))

"Usage: 
runDada2.R [options]

Description:   Run DADA2 on a Forward/Reverse fastq pair
Options:
--forwardReadpath=<forwardpath>                 Forward FastQ of paired end
--reverseReadpath=<reversepath>                 Forward FastQ of paired end
--samplename=<samplename>                       SampleName
--outdir=<outdir>               [default: '.']  Output Directory
--trimLeftf=<trimf>             [default: 0]    Number of bps to trim from forward read
--trimleftr=<trim>              [default: 0]    Number of bps to trim from reverse read
--truncLenf=<trunf>             [default: 225]  Number of bp at which to truncate the forward read
--truncLenr=<trunr>             [default: 150]  Number of bp at which to truncate the reverse read
--savetrimmed=<savetrimmed>     [default: TRUE] Save the trimmed fastqfiles?
--saveDADA2=<savadada>          [default: TRUE] Save the DADA2 fastqfiles?
--minfilesize=<minfilesize>     [default: 1000] Minimum filesize of DADA2 output.
--justConcatenate=<justconcat>  [default: TRUE] Concatenate? Merge is the default
--merge=<merge>                 [default: TRUE] Merge or Concatenate the reads
" -> doc
opts <- docopt(doc)

# check opts
fqf             <- opts$forwardReadpath
fqr             <- opts$reverseReadpath
samplename      <- opts$samplename
outdir          <- opts$outdir
trimLeftf       <- as.numeric(opts$trimLeftf)
trimLeftr       <- as.numeric(opts$trimLeftr)
truncLenf       <- as.numeric(opts$truncLenf)
truncLenr       <- as.numeric(opts$truncLenr)
minfilesize     <- as.numeric(opts$minfilesize)
savetrimmed     <- ifelse(opts$savetrimmed == "TRUE", TRUE, FALSE)
saveDADA2       <- ifelse(opts$saveDADA2 == "TRUE", TRUE, FALSE)
justConcatenate <- ifelse(opts$savetrimmed == "TRUE", TRUE, FALSE)
merge           <- ifelse(opts$saveDADA2 == "TRUE", TRUE, FALSE)

print(fqf)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dada2))


readPathsToDada2 <- function(forwardReadpath, 
                             reverseReadpath, 
                             samplename, 
                             outdir,
                             trimLeftf = 0,
                             trimLeftr = 0,
                             truncLenf = 225,
                             truncLenr = 150,
                             savetrimmed = TRUE,
                             saveDADA2 = TRUE,
                             minfilesize = 1000,
                             justConcatenate = TRUE,
                             merge=TRUE) {
  
  # deal with unpredicatable multithreading
  nthreads <- .Call(ShortRead:::.set_omp_threads, 1L)
  on.exit(.Call(ShortRead:::.set_omp_threads, nthreads))
  
  fqf_trimmed <- paste0(outdir,"/",samplename, "_F_filt_", truncLenf, ".fastq") 
  fqr_trimmed <- paste0(outdir,"/",samplename, "_R_filt_", truncLenr, ".fastq") 
  print(paste("Tempfiles: ", fqf_trimmed, fqr_trimmed))
  
  #' trim files
  filtered <- try(fastqPairedFilter(c(forwardReadpath, reverseReadpath), 
                                    c(fqf_trimmed,fqr_trimmed),
                                    maxN=0, 
                                    maxEE=2, 
                                    truncQ=2, 
                                    trimLeft=c(trimLeftf, trimLeftr), 
                                    truncLen=c(truncLenf,truncLenr), 
                                    compress=TRUE, 
                                    verbose=TRUE))
  print(file.info(fqf_trimmed)$size)
  print(file.info(fqf_trimmed)$size)
  
  #check trimmed files for size. If either is small, remove them
  if (file.info(fqf_trimmed)$size < minfilesize | file.info(fqr_trimmed)$size < minfilesize) {
    print(paste0("Outputfiles too small for sample ", samplename))
    file.remove(fqf_trimmed)
    file.remove(fqr_trimmed)
    return(NULL)
  }
  
  # dereplicate, remove file, perform dada, chimera checks, and merge 
  derepF <- derepFastq(fqf_trimmed, verbose=TRUE)
  derepR <- derepFastq(fqr_trimmed, verbose=TRUE)
  
  if (!savetrimmed) {
    file.remove(fqf_trimmed)
    file.remove(fqr_trimmed)
  }
  
  dadaF <- dada(derepF, err=inflateErr(tperr1,3), errorEstimationFunction=loessErrfun, selfConsist = TRUE)
  dadaR <- dada(derepR, err=inflateErr(tperr1,3), errorEstimationFunction=loessErrfun, selfConsist = TRUE)
  bimFs <- isBimeraDenovo(dadaF)
  bimRs <- isBimeraDenovo(dadaR)
  
  if (saveDADA2) {
    saveRDS(dadaF, file = paste0(outdir,"/",samplename,"_dadaF.RDS"))
    saveRDS(dadaR, file = paste0(outdir,"/",samplename,"_dadaR.RDS"))
    saveRDS(bimFs, file = paste0(outdir,"/",samplename,"_dadaF_bimera.RDS"))
    saveRDS(bimRs, file = paste0(outdir,"/",samplename,"_dadaR_bimera.RDS"))
    print(paste0("Saved DADA2 and bimera information for sample ",samplename))
  }
  
  if (merge) {
    mergers <- mergePairs(dadaF, derepF, dadaR, derepR, justConcatenate=justConcatenate)
    mergers.nochim <- mergers[!bimFs[mergers$forward] & !bimRs[mergers$reverse],]
    saveRDS(mergers.nochim,  file = paste0(outdir,"/",samplename,"_merged_", truncLenf, "_",truncLenr,".RDS"))
  }
}

readPathsToDada2(forwardReadpath = fqf,
                 reverseReadpath = fqr,
                 samplename = samplename, 
                 outdir = outdir,
                 trimLeftf = trimLeftf,
                 trimLeftr = trimLeftr,
                 truncLenf = truncLenf,
                 truncLenr = truncLenr,
                 savetrimmed = savetrimmed,
                 saveDADA2 = saveDADA2,
                 minfilesize = minfilesize,
                 justConcatenate = justConcatenate,
                 merge = merge)