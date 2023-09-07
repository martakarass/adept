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
  ind_mat[ind_mat > n] = NA_integer_
  ind_mat
}

make_shift_matrix = function(temp, nc) {
  stopifnot(abs(mean(temp)) <= 1e-5,
            abs(sd(temp) - 1) <= 1e-5)
  template_length = length(temp)
  last_value = temp[template_length]
  # zero pad it
  temp = c(temp, rep(0, nc - template_length))
  system.time({
    # quicker way for Toeplitz
    tt = c(temp[1], rep(0, nc - 1))
    shift_mat = pracma::Toeplitz(a = temp, b = tt)
    shift_mat[, (nc - template_length + 2):ncol(shift_mat)] = 0
  })
  stopifnot(
    shift_mat[nrow(shift_mat), (nc - template_length + 1)] == last_value
  )
  # shift_mat = shift_mat[, 1:(nc - template_length + 1)]
  shift_mat = Matrix::Matrix(shift_mat, sparse = TRUE)
}


make_shift_ones = function(template_length, nc) {
  temp = rep(TRUE, template_length)
  make_shift_matrix(temp, nc)
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

  x.smoothed = round(x.smoothed, 6)

  ind_mat = make_index_mat(x.smoothed, x.cut.vl, x.cut.seq, x.cut.margin)
  x_mat = array(x.smoothed[ind_mat], dim = dim(ind_mat))
  rm(ind_mat)
  nc = ncol(x_mat)
  not_na_x = !is.na(x_mat)
  not_na_x = Matrix::Matrix(not_na_x, sparse = TRUE)
  x_mat = Matrix::Matrix(x_mat, sparse = TRUE)



  template_list = template.scaled[[1]]

  result = pbmcapply::pbmclapply(template.scaled, function(template_list) {
    # result = pbapply::pblapply(template.scaled, function(template_list) {
    lengths = sapply(template_list, length)
    stopifnot(all(lengths == lengths[1]))
    template_length = lengths[1]

    one_mat = make_shift_ones(template_length, nc)
    n_mat = (not_na_x) %*% one_mat
    n_mat[n_mat <= 1L] = NA_integer_

    if (similarity.measure == "cor") {
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
    # result = pbapply::pblapply(template_list, function(temp) {
    res = lapply(template_list, function(temp) {
      shift_mat = make_shift_matrix(temp, nc)
      sum_xy = x_mat %*% shift_mat

      measure = sum_xy / denominator
      measure[!is.finite(measure)] = NA
      measure
    })
    res$na.rm = TRUE
    res = do.call(pmax, res)
    res
  }, mc.cores = getOption("mc.cores", mc.cores.val))

  reshaped = vector(mode = "list", length = nrow(x_mat))
  for (ireshaped in seq_along(reshaped)) {
    mat = matrix(nrow = length(result), ncol = nc)
    for (i in seq_along(result)) {
      mat[i, ] = result[[i]][ireshaped, ]
    }
    reshaped[[ireshaped]] = mat
  }
  rm(result)


  # back to x, not x.smoothed
  x = round(x, 6)
  finetune.maxima.x = round(finetune.maxima.x, 6)
  ind_mat = make_index_mat(x, x.cut.vl, x.cut.seq, x.cut.margin)
  x_mat = array(x[ind_mat], dim = dim(ind_mat))
  finetune_x_mat = array(finetune.maxima.x[ind_mat],
                         dim = dim(ind_mat))
  rm(ind_mat)

  final_results = pbmcapply::pbmclapply(
    seq_along(reshaped),
    function(i) {
      # for (i in seq_along(reshaped)) {
      # print(i)
      ii = x.cut.seq[i]
      xvals = x_mat[i, ]
      ft_vals = finetune_x_mat[i, ]
      similarity.mat = reshaped[[i]]
      if (i == length(reshaped)) {
        # end needs to be truncated
        na_x = is.na(xvals)
        xvals = xvals[!na_x]
        ft_vals = ft_vals[!na_x]
        similarity.mat = similarity.mat[,!na_x]
      }
      similarity.mat = round(similarity.mat, 8)
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
      ## Shift \tau parameter according to which part of signal x we are currently working with

      if (i == 1 || i == length(reshaped)) {
        ## Define current x part indices
        idx.i <- ii : min((ii + x.cut.vl + x.cut.margin), length(x.smoothed))
        ## If we cannot fit the longest pattern, return NULL
        if (length(idx.i) <= max(template.vl)) stop("FAIL")
        ## Compute similarity matrix
        similarity.mat.i <- similarityMatrix(
          x = round(x.smoothed[idx.i], 5),
          template.scaled = template.scaled,
          similarity.measure = similarity.measure)
        similarity.mat.i = round(similarity.mat.i, 8)
        stopifnot(isTRUE(
          all.equal(similarity.mat, similarity.mat.i, tolerance = 1e-5))
        )
        template.idx.mat.i <- NULL
        ## Run max and tine procedure
        out.df.i.check <- adept:::maxAndTune(
          x = round(x[idx.i], 5),
          template.vl = template.vl,
          similarity.mat = similarity.mat.i,
          similarity.measure.thresh = similarity.measure.thresh,
          template.idx.mat = template.idx.mat.i,
          finetune = finetune,
          finetune.maxima.x = round(finetune.maxima.x[idx.i], 5),
          finetune.maxima.nbh.vl = finetune.maxima.nbh.vl)
        stopifnot(isTRUE(all.equal(out.df.i, out.df.i.check, tolerance = 1e-5)))
      }

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
