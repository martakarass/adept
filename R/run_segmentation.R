# code taken from segmentPattern (old code/gold standard)
run_segmentation = function(x,
                            x.smoothed,
                            finetune.maxima.x,
                            x.cut.seq,
                            x.cut.vl,
                            x.cut.margin,
                            template.vl,
                            template.scaled,
                            similarity.measure,
                            similarity.measure.thresh,
                            compute.template.idx,
                            finetune,
                            finetune.maxima.nbh.vl,
                            mc.cores.val
) {

  out.list <- pbmcapply::pbmclapply(x.cut.seq, function(i){
    ## Define current x part indices
    idx.i <- i : min((i + x.cut.vl + x.cut.margin), length(x))
    ## If we cannot fit the longest pattern, return NULL
    if (length(idx.i) <= max(template.vl)) return(NULL)
    ## Compute similarity matrix
    similarity.mat.i <- similarityMatrix(x = x.smoothed[idx.i],
                                         template.scaled = template.scaled,
                                         similarity.measure = similarity.measure)
    ## Compute template index matrix
    if (compute.template.idx){
      template.idx.mat.i <- templateIdxMatrix(x = x.smoothed[idx.i],
                                              template.scaled = template.scaled,
                                              similarity.measure = similarity.measure)
    }
    else
    {
      template.idx.mat.i <- NULL
    }
    ## Run max and tine procedure
    out.df.i <- maxAndTune(x = x[idx.i],
                           template.vl = template.vl,
                           similarity.mat = similarity.mat.i,
                           similarity.measure.thresh = similarity.measure.thresh,
                           template.idx.mat = template.idx.mat.i,
                           finetune = finetune,
                           finetune.maxima.x = finetune.maxima.x[idx.i],
                           finetune.maxima.nbh.vl = finetune.maxima.nbh.vl)

    ## Shift \tau parameter according to which part of signal x we are currently working with
    if (nrow(out.df.i) > 0){
      out.df.i$tau_i <- out.df.i$tau_i + i - 1
      return(out.df.i)

    } else {
      ## Return empty data frame
      return(data.frame(tau_i = numeric(),
                        T_i = numeric(),
                        sim_i = numeric(),
                        template_i = numeric()))
    }
  }, mc.cores = getOption("mc.cores", mc.cores.val))
}
