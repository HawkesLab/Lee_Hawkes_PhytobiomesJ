# custom functions used in the dada2 pipeline

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found given a fastq path
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

all.primerHits.sample1 <- function(FWD, REV, Fs, Rs, sample.num){

  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)
    
  primerHitsTab <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = Fs[[sample.num]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = Rs[[sample.num]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = Fs[[sample.num]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = Rs[[sample.num]]))
  
  return(primerHitsTab)
}

all.adapterHits.sample1 <- function(adapter, Fs, Rs, sample.num){
  
  adapt.orients <- allOrients(adapter)
  hits.ForwardReads <- sapply(adapt.orients, primerHits, fn = Fs[[sample.num]])
  hits.ReverseReads <- sapply(adapt.orients, primerHits, fn = Rs[[sample.num]])
  result <- c(hits.ForwardReads, hits.ReverseReads)
  
  return(result)
}

readsLost_sample1 <- function(fastq.in, fastq.out, sample.num){
  
  pre.reads <- length(readFastq(fastq.in[sample.num]))
  post.reads <- length(readFastq(fastq.out[sample.num]))
  diff.reads <- pre.reads - post.reads
  lost.perc <- (diff/pre.reads) * 100
  
  data.frame(pre.reads = pre.reads, 
       post.reads = post.reads, 
       diff.reads = diff.reads, 
       lost.perc = lost.perc)
  
}

getN <- function(x) sum(getUniques(x))

get_nReads.seq <- function(x){
  nReads = as.numeric(getUniques(x))
  df <- data.frame(nReads =  nReads,
                   forward =seq(1:length(nReads)))
  return(df)
}

trackReads_old <- function(out, dadaFs, dadaRs, mergers, seqtab.nochim){
  
  ## SANITY CHECK: Track reads through the pipeline**********
  ## Look at the number of reads that made it through each step in the pipeline:
  ## Produce table where rows are samples, cols are num seqs at each step in processing
  ## Common to lose 5-10% of reads at each stage, but should never lose more than 40%
  
  #notes
  ## - If a majority of reads failed to merge, you may need to revisit `truncLen` to ensure overlap.
  ## - If a majority of reads were removed as chimeric, you may have unremoved primers.
  ## **This is the single most important place to inspect your workflow to make sure everything went as expected!**
  
  
  
  # before and after reads have been trimmed and filtered
  data.frame(out, row.names = NULL) %>%
    dplyr::rename("input"="reads.in",
                  "filtered"="reads.out") -> trim.filt
  
  # denoised, merged, chimeras removed
  denoised.merged.nonchim <- data.frame(
    denoisedF = sapply(dadaFs, getN),
    denoisedR = sapply(dadaRs, getN), 
    merged = sapply(mergers, getN),
    nonchim = rowSums(seqtab.nochim),
    row.names = NULL)
  
  tallytab <- data.frame(
    sample.name.match = row.names(seqtab.nochim),
    trim.filt,
    denoised.merged.nonchim)
  
  tallytab %>%
    mutate(perclost.filtered = (input-filtered)/input * 100) %>%
    mutate(perclost.denoisedF = (filtered-denoisedF)/filtered * 100) %>%
    mutate(perclost.denoisedR = (filtered-denoisedR)/filtered * 100) %>%
    mutate(perclost.mergedF = (denoisedF-merged)/denoisedF * 100) %>%
    mutate(perclost.mergedR = (denoisedR-merged)/denoisedR * 100) %>%
    mutate(perclost.nonchim = (merged-nonchim)/merged * 100) -> tallytab
  
  # plot
  tallytab %>%
    select(sample.name.match, input, filtered, denoisedF, denoisedR, merged, nonchim) %>%
    gather(key = "step",value = "numReads", -sample.name.match) -> tmp
  p <- ggplot(tmp, aes(x = sample.name.match, y = numReads, 
                       color = step)) +
    geom_point() +
    coord_flip() +
    xlab("Sample") + ylab("Number of reads")
  
  result <- list(tab = tallytab, p = p)
  
  return(result)
}

trackReads <- function(out, mergers, seqtab.nochim){
  
  ## SANITY CHECK: Track reads through the pipeline**********
  ## Look at the number of reads that made it through each step in the pipeline:
  ## Produce table where rows are samples, cols are num seqs at each step in processing
  ## Common to lose 5-10% of reads at each stage, but should never lose more than 40%
  
  #notes
  ## - If a majority of reads failed to merge, you may need to revisit `truncLen` to ensure overlap.
  ## - If a majority of reads were removed as chimeric, you may have unremoved primers.
  ## **This is the single most important place to inspect your workflow to make sure everything went as expected!**
  
  # before and after reads have been trimmed and filtered
  data.frame(out, row.names = NULL) %>%
    dplyr::rename("input"="reads.in",
                  "filtered"="reads.out") -> trim.filt
  
  # denoised, merged, chimeras removed
  merged.nonchim <- data.frame(
    merged = sapply(mergers, getN),
    nonchim = rowSums(seqtab.nochim),
    row.names = NULL)
  
  tallytab <- data.frame(
    sample.name.match = row.names(seqtab.nochim),
    trim.filt,
    merged.nonchim)
  
  tallytab %>%
    mutate(perclost.filtered = (input-filtered)/input * 100) %>%
    mutate(perclost.mergedF = (filtered-merged)/filtered * 100) %>%
    mutate(perclost.nonchim = (merged-nonchim)/merged * 100) -> tallytab
  
  return(tallytab)
}

plot_trackReads <- function(tracked){
  
  tracked %>%
    select(sample.name.match,
           perclost.filtered, perclost.mergedF, perclost.nonchim) %>%
    gather(key = "step", value = "perc", -sample.name.match) %>%
    separate(step, into = c("drop","step")) %>%
    select(-drop) %>%
    separate(sample.name.match, into = c("tissue","drop"), sep = 1, remove = F) %>%
    select(-drop) -> tmp
  step.levels <- c("filtered","mergedF","nonchim")
  tmp$step <- factor(tmp$step, levels = step.levels)
  tmp %>%
    filter(!grepl("NEG", sample.name.match)) -> tmp.nonegs
  
  p <- ggplot(tmp.nonegs, aes(x = step, y = perc, group = sample.name.match)) +
    geom_point(alpha = .5, pch = 16) +
    ylab("Percent of reads lost after each step") +
    facet_wrap(~tissue) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(p)
}

plot_initFinalReads <- function(tracked){
  
  tracked %>%
    select(sample.name.match,
           input, nonchim) %>%
    gather(key = "step", value = "reads", -sample.name.match) %>%
    separate(sample.name.match, into = c("tissue","drop"), sep = 1, remove = F) %>%
    select(-drop) -> tmp
  
  tmp %>%
    group_by(tissue, step) %>%
    summarize(n = length(sample.name.match),
              median = median(reads),
              lower25 = quantile(reads)[2],
              upper75 = quantile(reads)[4]) -> summ
  
  summ %>%
    filter(tissue %in% c("L","R","S")) -> summ.p
  p <- ggplot(summ.p, aes(x = tissue, y = median, color = step)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower25, ymax = upper75)) +
    geom_hline(yintercept = 10000, linetype = 2) +
    geom_hline(yintercept = 20000, linetype = 2) +
    xlab("Sample compartment") + ylab("Reads per sample (lower 25%, median, upper 75%)") +
    theme_bw()
  
  result <- list(p = p, summ = summ)
  return(result)
}

examinePairedReads <- function(mergers.obj, dadaFs.obj, dadaRs.obj, mergers.summ, 
                               curr.sample, F.indx, R.indx){
  
  # forward read
  F.read <- DNAStringSet(names(dadaFs.obj[[curr.sample]]$denoised[F.indx]))
  
  # reverse read (in reverse complement format)
  R.read <- DNAStringSet(names(dadaRs.obj[[curr.sample]]$denoised[R.indx]))
  R.read.rc <- reverseComplement(R.read)
  
  # merged sequence
  selection <- mergers.summ$sample.name.match == curr.sample &
    mergers.summ$forward == F.indx &
    mergers.summ$reverse == R.indx
  summ.row <- which(selection)
  merged.seq <- DNAStringSet(mergers.obj[[curr.sample]]$sequence[[summ.row]])
  
  
  seqs <- c(F.read, R.read.rc, merged.seq)
  names(seqs) <- c("forward","reverse.rc","merged")
  
  return(seqs)
  
}

addFR_lengths <- function(dadaFs.obj, dadaRs.obj, mergers.summ){
  
  # annotate table with length of forward and reverse reads
  
  # forward reads -- include the whole sequence too
  lengthFs.list <- list()
  for(i in 1:length(dadaFs.obj)){
    length.f <- nchar(names(dadaFs.obj[[i]]$denoised))
    lengthFs.list[[i]] <- data.frame(forward = seq(1, length(length.f)), 
                                     length.f = length.f,
                                     seq.f = names(dadaFs.obj[[i]]$denoised))
  }
  names(lengthFs.list) <- names(dadaFs.obj)
  lengthFs.df <- list_to_df(lengthFs.list)
  lengthFs.df %>%
    dplyr::rename("sample.name.match"="source") -> lengthFs.df
  
  # reverse reads
  lengthRs.list <- list()
  for(i in 1:length(dadaFs.obj)){
    length.r <- nchar(names(dadaRs.obj[[i]]$denoised))
    lengthRs.list[[i]] <- data.frame(reverse = seq(1, length(length.r)), 
                                     length.r = length.r)
  }
  names(lengthRs.list) <- names(dadaRs.obj)
  lengthRs.df <- list_to_df(lengthRs.list)
  lengthRs.df %>%
    dplyr::rename("sample.name.match"="source") -> lengthRs.df
  
  #annotate df
  mergers.summ %>%
    left_join(lengthFs.df) %>%
    left_join(lengthRs.df) %>%
    mutate(length.sum = length.f + length.r) -> mergers.summ.ann
  
  return(mergers.summ.ann)
}

make_f.denoised_df <- function(dadaFs.obj){
  # all forward seqs
  seqFs.list <- list()
  for(i in 1:length(dadaFs.obj)){
    seq.f <- names(dadaFs.obj[[i]]$denoised)
    length.f <- nchar(seq.f)
    abundance.f <- getUniques(dadaFs.obj[[i]]$denoised)
    seqFs.list[[i]] <- data.frame(forward = seq(1, length(length.f)), 
                                  length.f = length.f,
                                  abundance.f = abundance.f,
                                  seq.f = seq.f)
  }
  names(seqFs.list) <- names(dadaFs.obj)
  seqFs.df <- list_to_df(seqFs.list)
  seqFs.df %>%
    dplyr::rename("sample.name.match"="source") -> seqFs.df
  
  return(seqFs.df)
}

make_unmerged.f_df <- function(mergers.summ, seqFs.df){
  
  mergers.summ %>%
    group_by(sample.name.match, forward, accept) %>%
    summarize(n = length(forward)) %>%
    spread(key = "accept",value = "n", fill = 0) %>%
    filter(`TRUE`==0) -> lost.f.seqs
  
  # add sequences to lost.f.seqs and add unique names
  lost.f.seqs %>%
    left_join(seqFs.df) %>%
    mutate(seq.name = paste(sample.name.match, forward, sep = "_")) -> lost.f.seqs
  
  return(lost.f.seqs)
}

make_mergeforward.df <- function(mergers.summ, dadaFs){
  
  # subset forward unmerged sequences per sample
  mergers.summ %>%
    filter(accept == FALSE) %>%
    select(forward, sample.name.match, length.f, seq.f) %>%
    unique() -> f.unmerged.df 
  
  # subset forward merged sequences per sample
  mergers.summ %>%
    filter(accept == TRUE) %>%
    select(forward, sample.name.match, length.f) %>%
    unique() %>%
    mutate(mergestat = "merged") -> f.merged.df
  
  # identify forwards that are not represented in at least 1 merged sequence
  f.unmerged.df %>%
    left_join(f.merged.df) %>%
    mutate(mergestat = ifelse(is.na(mergestat),
                              "unmerged", mergestat)) -> f.seqs.df
  
  # annotate with abundance from previous step
  dadaFs %>%
    map(~get_nReads.seq(.x)) %>%
    list_to_df %>%
    dplyr::rename('sample.name.match'='source') -> denoisedF.df
  f.seqs.df %>%
    left_join(denoisedF.df) -> f.seqs.df.ann
  
  return(f.seqs.df.ann)
}

perc_unmergedforward <- function(dadaFs, mergeforward){
  
  # summarize denoised forwards
  dadaFs %>%
    map(~getN(.x)) -> nreads.denoisedf
  nreads.denoisedf.df <- data.frame(sample.name.match = names(nreads.denoisedf),
                                    denoisedf = as.numeric(nreads.denoisedf))
  
  # summarize the length of forward merged and unmerged sequences per sample, add denoised forwards
  mergeforward %>%
    group_by(sample.name.match, mergestat) %>%
    summarize(nReads = sum(nReads)) %>%
    spread(key = "mergestat",value = "nReads") %>%
    left_join(nreads.denoisedf.df) -> df
  
  # summarize the percent lost to merge
  df %>%
    mutate(perc.unmergedOfDenoised = (unmerged/denoisedf)*100) %>%
    mutate(perc.mergedOfDenoised = (merged/denoisedf)*100) -> df.summ
  
  return(df.summ)
  
}

summarize_read.length <- function(mergeforward, mergers.summ){
  
  # summarize the forward reads
  mergeforward %>%
    group_by(sample.name.match, mergestat) %>%
    summarize(n = length(length.f),
              minlength = min(length.f),
              maxlength = max(length.f),
              medianlength = median(length.f)) %>%
    mutate(read.type = ifelse(mergestat == "merged",
                              "merged_forward",
                              "unmerged_forward")) %>%
    select(-mergestat) -> forwards.df
  
  # summarize the merged reads
  mergers.summ %>%
    filter(accept == TRUE) %>%
    group_by(sample.name.match) %>%
    summarize(n = length(nchar.merge.seq),
              minlength = min(nchar.merge.seq),
              maxlength = max(nchar.merge.seq),
              medianlength = median(nchar.merge.seq)) %>%
    mutate(read.type = "merged") -> merged.df
  
  # combine
  forwards.df %>%
    bind_rows(merged.df) %>%
    arrange(sample.name.match, read.type) %>%
    select(read.type, sample.name.match, n, minlength, maxlength, medianlength) -> length.df
  
  return(length.df)
}

pull_long.unmergedfs <- function(lengths, mergeforward, add.bases){
  # find the median max length of merged forwards
  lengths %>%
    filter(read.type == "merged_forward") -> tmp
  max.len <- median(tmp$maxlength)
  
  # subset unmerged forwards that are at least 3bps longer than this
  mergeforward %>%
    filter(mergestat == "unmerged") %>%
    filter(length.f > max.len + add.bases) -> long.unmergedfs
  
  return(long.unmergedfs)
}
