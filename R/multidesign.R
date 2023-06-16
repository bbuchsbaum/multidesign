
#' @importFrom dplyr mutate rowwise n
#' @export
#' @examples
#'
#' X <- matrix(rnorm(20*100), 20, 100)
#' Y <- tibble(condition=rep(letters[1:5], 4), subject=rep(1:4, each=5))
#'
#' mds <- multidesign(X,Y)
#' @rdname multidesign
multidesign.matrix <- function(x, y) {
  chk::chk_equal(nrow(x), nrow(y))
  chk::chk_s3_class(y, "data.frame")
  y <- as_tibble(y)

  des <- y %>% mutate(.index=1:n())
  structure(list(
    x=x,
    design=tibble::as_tibble(des)
  ),
  class="multidesign")
}

reduce.multidesign <- function(x, nc=2, ..., rfun=function(x) multivarious::pca(x$x, ncomp=nc,...)) {
  projector <- rfun(x)
  chk::chk_s3_class(projector, "projector")
  rx <- multivarious::project(projector, x$x)
  structure(list(
    x=rx,
    design=x$design,
    projector=projector
  ),
  class=c("reduced_multidesign", "multidesign"))
}

#' @export
subset.multidesign <- function(x, fexpr) {
  des2 <- filter(x$design, !!rlang::enquo(fexpr))
  if (nrow(des2) == 0) {
    NULL
  } else {
    ind <- des2$.index
    multidesign(x$x[ind,], des2)
  }
}

#' @export
split.multidesign <- function(x, ..., collapse=FALSE) {
  nest.by <- rlang::quos(...)
  ret <- x$design %>% nest_by(!!!nest.by, .keep=TRUE)
  xl <- ret$data %>% purrr::map(~x$x[.x$.index,,drop=FALSE])
  ret <- lapply(1:nrow(ret), function(i) multidesign.matrix(xl[[i]], ret$data[[i]]))
}

#' @export
#' @importFrom tidyr unite
split_indices.multidesign <- function(x, ..., collapse=FALSE) {
  nest.by <- rlang::quos(...)
  ret <- x$design %>% nest_by(!!!nest.by, .keep=TRUE)
  xl <- ret$data %>% purrr::map(~ .x$.index)
  cvars <- ret %>% select(group_vars(ret)) %>% unite(nest.by) %>% pull("nest.by")
  selret <- ret %>% select(group_vars(ret))
  out <- selret %>% ungroup() %>% mutate(indices=xl, .splitvar=cvars)
  out
}


#' @export
summarize_by.multidesign <- function(x, ..., sfun=colMeans, extract_data=FALSE) {
  #nested <- split(x, ...)
  nest.by <- rlang::quos(...)
  ret <- x$design %>% nest_by(!!!nest.by)

  dsum <- do.call(rbind, ret$data %>% purrr::map( ~ sfun(x$x[.x[[".index"]], drop=FALSE,])))
  ret2 <- ret %>% select(-data)
  multidesign(dsum, ret2)
}

#' @export
xdata.multidesign <- function(x) x$x

#' @export
design.multidesign <- function(x) x$y


#' @export
#' @importFrom deflist deflist
fold_over.multidesign <- function(x, ...) {
  args <- rlang::enquos(...)
  splits <- split_indices(x,!!!args)
  foldframe <- splits %>% mutate(.fold=1:n())

  extract <- function(i) {
    #browser()
    block <- foldframe[[".fold"]][i]
    ind <- unlist(foldframe[["indices"]][[i]])

    testdat <- multidesign(xdata(x)[ind,], x$design[ind,])
    ## all blocks except
    traindat <- multidesign(xdata(x)[-ind,], x$design[-ind,])
    list(analysis=traindat,
         assessment=testdat)

  }

  tlen <- nrow(foldframe)
  ret <- deflist(extract, len=tlen)
  names(ret) <- paste0("fold_", 1:length(ret))
  class(ret) <- c("foldlist", class(ret))
  attr(ret, "foldframe") <- foldframe
  ret
}


#' @export
print.multidesign <- function(x) {
  cat("A multidesign object. \n")
  cat(nrow(x$design), "rows", "\n")
  cat(ncol(x$des), "design variables", "\n")
  cat("design variables: ", "\n")
  print(x$design, n=5)
}

#' @export
print.reduced_multidesign <- function(x) {
  cat("A (reduced) multidesign object. \n")
  cat(paste("original input channels: ", multivarious::shape(x$projector)[1], "\n"))
  cat(paste("reduced input channels: ", multivarious::shape(x$projector)[2], "\n"))
  cat(nrow(x$design), "rows", "\n")
  cat(ncol(x$design), "design variables")
  cat("design variables: ", "\n")
  print(x$design, n=5)

}
