
# Make an index matrix for the window
make_index_mat = function(data, x.cut.vl, x.cut.seq, x.cut.margin) {
  n = length(data)
  nr = length(x.cut.seq)
  ind_mat = matrix(1:(nr * x.cut.vl), ncol = x.cut.vl, nrow = nr,
                   byrow = TRUE)
  addon_mat = matrix(rep(1:(x.cut.margin + 1), times = nr),
                     nrow = nr,
                     ncol = x.cut.margin + 1,
                     byrow = TRUE)
  addon_mat = addon_mat + ind_mat[, ncol(ind_mat)]
  ind_mat = cbind(ind_mat, addon_mat)

  # make NA indices that are longer than the vecor
  ind_mat[ind_mat > n] = NA_integer_
  ind_mat
}

make_shift_index_matrix = function(nc) {
  suppressWarnings({
    res = matrix(nrow = nc, ncol = nc)
    res[1:(nc-1), 1:(nc-1)] = matrix(1:(nc), nc-1, nc-1)
    res[, nc] = c(2:nc, 1)
    res[nc, ] = nc:1
  })
  res
}

# create the rolling shift matrix
# should be
# x1 0 0 ...
# x2 x1 0 ...
# x3 x2 x1 0 ...
make_shift_matrix = function(temp, nc, check_values = TRUE ) {
  if (check_values) {
    stopifnot(abs(mean(temp)) <= 1e-5,
              abs(stats::sd(temp) - 1) <= 1e-5)
  }
  template_length = length(temp)
  last_value = temp[template_length]
  # zero pad it
  temp = c(temp, rep(0, nc - template_length))
  ind_mat = make_shift_index_matrix(nc = nc)
  shift_mat = array(temp[ind_mat], dim = dim(ind_mat))
  # quick check
  stopifnot(
    shift_mat[nrow(shift_mat), (nc - template_length + 1)] == last_value
  )
  shift_mat = Matrix::Matrix(shift_mat, sparse = TRUE)
}


# creates the same thing shifted, but just TRUE (1s) and FALSE (0)
# this is to get rolling sums for X and X^2
make_shift_ones = function(template_length, nc) {
  temp = rep(TRUE, template_length)
  make_shift_matrix(temp, nc, check_values = FALSE)
}


run_fast_segmentation = function(
    x,
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

  # round the data for numerical stability
  x.smoothed = round(x.smoothed, 6)

  ind_mat = make_index_mat(x.smoothed, x.cut.vl, x.cut.seq, x.cut.margin)
  # get all the sequences out
  x_mat = array(x.smoothed[ind_mat], dim = dim(ind_mat))
  rm(ind_mat)
  nc = ncol(x_mat)
  # will use this to get n for the cov(x,y)/n and cor(x, y)
  not_na_x = !is.na(x_mat)
  not_na_x = Matrix::Matrix(not_na_x, sparse = TRUE)
  x_mat = Matrix::Matrix(x_mat, sparse = TRUE)



  template_list = template.scaled[[1]]

  # using pbmclapply for a progress bar
  result = pbmcapply::pbmclapply(template.scaled, function(template_list) {
    # result = pbapply::pblapply(template.scaled, function(template_list) {

    # quick template checks.  These are crucial because one speedup
    # is that we only need to calculate the values for x once
    # per template stretch list
    lengths = sapply(template_list, length)
    stopifnot(all(lengths == lengths[1]))
    template_length = lengths[1]

    # get n for the cov/cor calculation
    index_mat = make_shift_index_matrix(nc)
    one_mat = c(rep(TRUE, template_length),
                rep(FALSE, nc - template_length))
    one_mat = array(one_mat[index_mat], dim = dim(index_mat))
    # one_mat = make_shift_ones(template_length, nc)
    n_mat = (not_na_x) %*% one_mat
    n_mat[n_mat <= 1L] = NA_integer_

    if (similarity.measure == "cor") {
      # matrix version of correlation, but rolling
      # has BIG assumption of E[y] = 0 and sd(y) = 1

      sum_x2 = (x_mat ^ 2) %*% one_mat
      sum_x = x_mat %*% one_mat
      # 1/(n-1) (Σ(x - 𝑥̄)^2
      # 1/(n-1) (Σx^2 - n𝑥̄^2)
      # 1/(n-1) (Σx^2 - n (Σx/n)^2)
      # 1/(n-1) (Σx^2 - (Σx)^2/n)
      # denominator = sqrt(1/(n_mat - 1) * (sum_x2 - (sum_x ^ 2)/n_mat))
      denominator = round(sum_x2 - (sum_x ^ 2)/n_mat, 8)
      denominator = sqrt(denominator)
      # { 1/(n-1) Σ(x_i y_i) } / √{(1/n-1) SS_x}
      # {  Σ(x_i y_i) } / √(n-1) √{SS_x}
      rm(sum_x)
      rm(sum_x2)
      denominator = denominator * sqrt(n_mat - 1L)
    } else {
      denominator = n_mat - 1L
    }
    rm(one_mat)

    # itemp = 1
    # temp = template_list[[1]]
    res = lapply(template_list, function(temp) {

      # here is the rolling cross product
      # shift_mat = make_shift_matrix(temp, nc)
      temp_pad = c(temp, rep(0, nc - length(temp)))
      shift_mat = array(temp_pad[index_mat], dim = dim(index_mat))

      # shift_mat = array(temp, dim = dim(index_mat))
      # shift_mat[is.na(shift_mat)] = 0
      sum_xy = x_mat %*% shift_mat

      measure = sum_xy / denominator
      measure[!is.finite(measure)] = NA
      measure
    })
    res$na.rm = TRUE
    # running pmax here.  We would need to do a workup on `res`
    # to get template.idx
    res = do.call(pmax, res)
    res
  }, mc.cores = getOption("mc.cores", mc.cores.val))

  # need to reshape it because maxAndTune requires a different list
  reshaped = vector(mode = "list", length = nrow(x_mat))
  for (ireshaped in seq_along(reshaped)) {
    mat = matrix(nrow = length(result), ncol = nc)
    for (i in seq_along(result)) {
      mat[i, ] = result[[i]][ireshaped, ]
    }
    reshaped[[ireshaped]] = mat
  }
  # remove the old result
  rm(result)


  # back to x, not x.smoothed
  # rounding for stability
  x = round(x, 6)


  # again making an index matrix
  ind_mat = make_index_mat(x, x.cut.vl, x.cut.seq, x.cut.margin)
  x_mat = array(x[ind_mat], dim = dim(ind_mat))
  run_finetune = !is.null(finetune) && finetune == "maxima"
  if (run_finetune) {
    finetune.maxima.x = round(finetune.maxima.x, 6)
    finetune_x_mat = array(finetune.maxima.x[ind_mat],
                           dim = dim(ind_mat))
  } else {
    finetune.maxima.x = NULL
    finetune_x_mat = NULL
  }

  # we could rm(x.smoothed) here but we use it for checking
  rm(ind_mat)

  final_results = pbmcapply::pbmclapply(
    seq_along(reshaped),
    function(i) {
      # for (i in seq_along(reshaped)) {
      # print(i)
      ii = x.cut.seq[i]
      xvals = x_mat[i, ]
      if (run_finetune) {
        ft_vals = finetune_x_mat[i, ]
      } else {
        ft_vals = NULL
      }
      similarity.mat = reshaped[[i]]

      # need workup for the last row because not always fit nicely
      # (e.g. end of vector)
      if (i == length(reshaped)) {
        # end needs to be truncated
        na_x = is.na(xvals)
        xvals = xvals[!na_x]
        if (run_finetune) {
          ft_vals = ft_vals[!na_x]
        }
        similarity.mat = similarity.mat[, !na_x]
      }
      # round this to something reasonable
      similarity.mat = round(similarity.mat, 8)

      # again, would need workup to make template.idx to work
      template.idx.mat.i = NULL
      out.df.i <- maxAndTune(
        x = xvals,
        template.vl = template.vl,
        similarity.mat = similarity.mat,
        similarity.measure.thresh = similarity.measure.thresh,
        template.idx.mat = template.idx.mat.i,
        finetune = finetune,
        finetune.maxima.x = ft_vals,
        finetune.maxima.nbh.vl = finetune.maxima.nbh.vl)

      # for first and last record, check against old version/gold standard
      if (i == 1 || i == length(reshaped)) {
        ## Define current x part indices
        idx.i <- ii : min((ii + x.cut.vl + x.cut.margin), length(x.smoothed))
        ## If we cannot fit the longest pattern, return NULL
        if (length(idx.i) <= max(template.vl)) stop("FAIL")
        ## Compute similarity matrix
        similarity.mat.i <- similarityMatrix(
          x = x.smoothed[idx.i],
          template.scaled = template.scaled,
          similarity.measure = similarity.measure)
        similarity.mat.i = round(similarity.mat.i, 8)
        # need some reasonable tolerance here
        stopifnot(isTRUE(
          all.equal(similarity.mat, similarity.mat.i, tolerance = 1e-5))
        )
        template.idx.mat.i <- NULL
        ## Run max and tine procedure
        out.df.i.check <- maxAndTune(
          x = x[idx.i],
          template.vl = template.vl,
          similarity.mat = similarity.mat.i,
          similarity.measure.thresh = similarity.measure.thresh,
          template.idx.mat = template.idx.mat.i,
          finetune = finetune,
          finetune.maxima.x = finetune.maxima.x[idx.i],
          finetune.maxima.nbh.vl = finetune.maxima.nbh.vl)
        stopifnot(isTRUE(all.equal(out.df.i, out.df.i.check, tolerance = 1e-5)))
      }

      # shift the tau based on ii - the cut index
      if (nrow(out.df.i) > 0){
        out.df.i$tau_i <- out.df.i$tau_i + ii - 1
      } else {
        ## Return empty data frame
        out.df.i = data.frame(tau_i = numeric(),
                              T_i = numeric(),
                              sim_i = numeric(),
                              template_i = numeric())
      }
      return(out.df.i)
    }, mc.cores = getOption("mc.cores", mc.cores.val))


  return(final_results)
}
