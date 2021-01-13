#' @keywords internal
block_index_mat <- function(x, byrow=FALSE) {
  lens <- if (byrow) {
    sapply(x, function(z) nrow(z$x))
  } else {
    sapply(x, function(z) ncol(z$x))
  }

  csum <- cumsum(lens)
  csum1 <- c(0, csum[-length(csum)])
  m <- as.matrix(cbind(csum1+1, csum))
  colnames(m) <- c("start", "end")
  m
}

#' @export
#'
#' @param block_names the names of each block
#'
#' @examples
#' d1 <- multidesign(matrix(rnorm(10*20), 10, 20), data.frame(y=1:10, subject=1, run=rep(1:5, 2)))
#' d2 <- multidesign(matrix(rnorm(10*20), 10, 20), data.frame(y=1:10, subject=2, run=rep(1:5, 2)))
#' d3 <- multidesign(matrix(rnorm(10*20), 10, 20), data.frame(y=1:10, subject=3, run=rep(1:5, 2)))
#'
#' hd <- hyperdesign(list(d1,d2,d3))
hyperdesign <- function(x, block_names=NULL) {
  chk::chk_true(all(sapply(x, function(d) inherits(d, "multidesign"))))

  bind_col <- block_index_mat(x, byrow=FALSE)
  bind_row <- block_index_mat(x, byrow=TRUE)

  if (!is.null(block_names)) {
    chk::chk_true(length(block_names) == length(x))
    names(x) <- block_names
  } else if (is.null(names(x))) {
    names(x) <- paste0("block_", 1:length(x))
  }

  hdes <- lapply(1:length(x), function(i) {
    tibble(block=i, nr=nrow(x[[i]]$x), nxvar=ncol(x[[i]]$x), nyvar=ncol(x[[i]]$design),
           row_start=bind_row[i,1], row_end=bind_row[i,2],
           col_start=bind_col[i,1], col_end=bind_col[i,2])
  }) %>% bind_rows()

  cvars <- lapply(x, function(z) names(z$design))
  cvars <- Reduce(intersect, cvars)

  structure(x,
            hdes=hdes,
            common_vars=cvars,
            class="hyperdesign")
}

#' @export
#' @param i the block number
#' @param byrow if true, return row-oriented indices
#' @rdname block_indices
block_indices.hyperdesign <- function(x, i, byrow=FALSE) {
  hd <- attr(x, "hdes")

  if (missing(i)) {
    lapply(1:nrow(hd), function(j) {
      if (byrow) {
        seq(hd$row_start[i], hd$row_end[i])
      } else {
        seq(hd$col_start[i], hd$col_end[i])
      }
    })
  } else {
    chk::chk_range(i, c(1, nrow(hd)))
    if (byrow) {
      seq(hd$row_start[i], hd$row_end[i])
    } else {
      seq(hd$col_start[i], hd$col_end[i])
    }
  }
}

#' @importFrom multivarious init_transform
#' @return a transformed `hyperdesign` object, with pre-processors set as attribute named `preproc`
#' @export
#' @rdname init_transform
init_transform.hyperdesign <- function(x, preproc) {
  ## pre-processors
  proclist <- lapply(seq_along(x), function(i) {
    multivarious:::fresh(preproc) %>% prep()
  })

  names(proclist) <- names(x)

  ## subject-split and pre-processed data
  out <- lapply(seq_along(proclist), function(i) {
    p <- proclist[[i]]
    Xi <- x[[i]]$x
    Xout <- multivarious::init_transform(p, Xi)
    multidesign(Xout, x[[i]]$design)
  })

  des <- hyperdesign(out)
  attr(des, "preproc") <- proclist
  des
}


#' @examples
#' d1 <- multidesign(matrix(rnorm(10*20), 10, 20), data.frame(y=1:10, subject=1, run=rep(1:5, 2)))
#' d2 <- multidesign(matrix(rnorm(10*20), 10, 20), data.frame(y=1:10, subject=2, run=rep(1:5, 2)))
#' d3 <- multidesign(matrix(rnorm(10*20), 10, 20), data.frame(y=1:10, subject=3, run=rep(1:5, 2)))
#'
#' hd <- hyperdesign(list(d1,d2,d3))
#' folds <- fold_over(hd, run)
#' @export
#' @rdname fold_over
fold_over.hyperdesign <- function(x, ...) {

  splits <- lapply(seq_along(x), function(i) {
    d <- x[[i]]
    split_indices(d, ...) %>% mutate(.block=i)
  })

  foldframe <- splits %>% bind_rows() %>% mutate(.fold=1:n())

  lens <- purrr::map(splits, function(x) nrow(x))

  tlen <- sum(unlist(lens))


  extract <- function(i) {
    block <- foldframe[[".block"]][i]
    ind <- unlist(foldframe[["indices"]][[i]])

    testdat <- multidesign(x[[block]]$x[ind,], x[[block]]$design[ind,])

    traindat <- hyperdesign(lapply(seq_along(x), function(j) {
      if (j == block) {
        if (length(ind) == nrow(x[[j]]$x)) {
          stop(paste("number of test data rows in block ", j, " is equal to number of rows in block.
                     No training data available for current fold! Check distribution of splitting variables over blocks."))
        }
        multidesign(x[[block]]$x[-ind,], x[[block]]$design[-ind,])
      } else{
        x[[block]]
      }
    }))

    list(analysis=traindat,
         assessment=testdat)

  }

  deflist(extract, len=tlen)
  ## split over the variable, get global indices
  ## verify that variable exists in each design
  ## verify that all blocks have a minimum observation size
  ## create a `deflist` ,wher each element has an "analysis" and assessment hyperdesign

}

#' @export
xdata.hyperdesign <- function(x, block) {
  if (missing(block)) {
    lapply(x, xdata)
  } else {
    chk::vld_number(block)
    chk::chk_range(block, 1, length(x))
    xdata(x[[block]])
  }
}

#' @export
design.hyperdesign <- function(x, block) {
  if (missing(block)) {
    lapply(x, design)
  } else {
    chk::vld_number(block)
    chk::chk_range(block, 1, length(x))
    design(x[[block]])
  }
}



#' @export
print.hyperdesign <- function(x) {
  cat("a hyperdesign object. \n")
  cat("with ", length(x), "blocks.", "\n")
  cat("common design variables: ", attr(x, "common_vars"), "\n")
  print(attr(x, "hdes"))
  cat("\n")

}
