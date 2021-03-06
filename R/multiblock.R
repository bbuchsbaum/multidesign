#' @keywords internal
block_index_mat2 <- function(x, byrow=FALSE) {
  lens <- if (byrow) {
    sapply(x, function(z) nrow(z))
  } else {
    sapply(x, function(z) ncol(z))
  }

  csum <- cumsum(lens)
  csum1 <- c(0, csum[-length(csum)])
  m <- as.matrix(cbind(csum1+1, csum))
  colnames(m) <- c("start", "end")
  m
}

#' @export
multiblock.list <- function(x, design=NULL) {
  sapply(x, function(m) chk::chkor(chk::chk_s3_class(m, "matrix"), chk::chk_s4_class(m, "Matrix")))

  nr <- sapply(x, nrow)
  nc <- sapply(x, ncol)

  ret <- if (all(nr == nr[1])) {
    ## c stacked
    bind <- block_index_mat2(x, byrow=TRUE)
    structure(x,
              ind=bind,
              orient="cstacked",
              class=c("multiblock_c", "multiblock_list", "multiblock", "list"))
  } else if (all(nc == nc[1])) {
    ## r stacked
    bind <- block_index_mat2(x, byrow=FALSE)
    structure(x,
              ind=bind,
              orient="rstacked",
              class=c("multiblock_r", "multiblock_list", "multiblock", "list"))
  }
}

#' @export
block_indices.multiblock_list <- function(x, i) {
  chk::chk_range(i, c(1,length(x)))
  ind <- attr(x, "ind")
  seq(ind[i, 1], ind[i,2])
}


#' @export
print.multiblock_list <- function(x) {
  cat("a multiblock object with", length(x), "blocks \n")
  #cat(length(x), "blocks", "\n")
  cat("orientation: ", attr(x, "orient"), "\n")
  if (attr(x, "orient") == "cstacked") {
    cat("rows: ", nrow(x[[1]]), "\n")
    cat("columns: ", paste(sapply(x, ncol), sep=" "))
  } else {
    cat("rows: ", paste(sapply(x,nrow), sep=" "), "\n")
    cat("columns: ", ncol(x[[1]]), "\n")
  }


}


