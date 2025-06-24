#' Create Index Matrix for Block Structure
#'
#' @description Internal function to compute start and end indices for blocks in a multiblock object
#' @param x List of matrices
#' @param byrow Logical; if TRUE compute row indices, if FALSE compute column indices
#' @return Matrix with two columns: 'start' and 'end' indices for each block
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

#' Create a Multiblock Object from a List of Matrices
#'
#' @description 
#' Creates a multiblock object from a list of matrices that share either row or column dimensions.
#' The function automatically determines whether the matrices should be row-stacked or column-stacked
#' based on their shared dimensions.
#'
#' @details
#' The function checks if the input matrices can be combined either by:
#' * Row-stacking: all matrices must have the same number of columns
#' * Column-stacking: all matrices must have the same number of rows
#'
#' The resulting object maintains block structure information while allowing operations
#' across the entire combined matrix.
#'
#' @param x A list of matrices (base or Matrix class objects)
#' @param ... Additional arguments (not used)
#' @return A multiblock object with the following attributes:
#'   * ind: matrix of start/end indices for each block
#'   * orient: orientation of stacking ("cstacked" or "rstacked")
#'   * class: appropriate class labels for dispatch
#'
#' @family multiblock functions
#' @seealso 
#'   [is_cstacked()] for checking stacking orientation,
#'   [block_indices()] for accessing block-specific indices
#'
#' @examples
#' # Create example matrices with shared row dimension (column-stacked)
#' X1 <- matrix(rnorm(20*3), 20, 3)
#' X2 <- matrix(rnorm(20*5), 20, 5)
#' X3 <- matrix(rnorm(20*4), 20, 4)
#' mb_c <- multiblock(list(X1, X2, X3))  # column-stacked
#' is_cstacked(mb_c)  # TRUE
#'
#' # Create example matrices with shared column dimension (row-stacked)
#' Y1 <- matrix(rnorm(10*5), 10, 5)
#' Y2 <- matrix(rnorm(15*5), 15, 5)
#' Y3 <- matrix(rnorm(20*5), 20, 5)
#' mb_r <- multiblock(list(Y1, Y2, Y3))  # row-stacked
#' is_rstacked(mb_r)  # TRUE
#' @export
multiblock.list <- function(x, ...) {
  if (length(x) == 0) {
    stop("List is empty")
  }
  
  # Check if all elements are matrices
  is_matrix <- sapply(x, function(m) inherits(m, "matrix") || inherits(m, "Matrix"))
  if (!all(is_matrix)) {
    stop("Not all elements are matrices")
  }

  nr <- sapply(x, nrow)
  nc <- sapply(x, ncol)

  ret <- if (all(nr == nr[1])) {
    ## c stacked
    bind <- block_index_mat2(x, byrow=FALSE)  
    structure(x,
              ind=bind,
              orient="cstacked",
              class=c("multiblock_c", "multiblock_list", "multiblock", "list"))
  } else if (all(nc == nc[1])) {
    ## r stacked
    bind <- block_index_mat2(x, byrow=TRUE)  
    structure(x,
              ind=bind,
              orient="rstacked",
              class=c("multiblock_r", "multiblock_list", "multiblock", "list"))
  } else {
    stop("all matrices must share either row or column dimension")
  }
}

#' Extract Block Indices from Multiblock Object
#'
#' @description
#' Retrieves the indices corresponding to a specific block in a multiblock object.
#' These indices can be used to access the corresponding rows or columns in the 
#' combined matrix representation.
#'
#' @param x A multiblock object
#' @param i Integer specifying which block's indices to retrieve
#' @param ... Additional arguments (not used)
#' @return Integer vector of indices for the specified block
#' @family multiblock functions
#' @method block_indices multiblock_list
#' @export
#' 
#' @examples 
#' mb <- multiblock(list(matrix(1:10, 5, 2), matrix(11:20, 5, 2)))
#' block_indices(mb, 1)
#' block_indices(mb, 2)
block_indices.multiblock_list <- function(x, i, ...) {
  chk::chk_range(i, c(1,length(x)))
  ind <- attr(x, "ind")
  seq(ind[i, 1], ind[i,2])
}

#' Transpose a Multiblock Object
#'
#' @description
#' Transposes each matrix in a multiblock object and returns a new multiblock object.
#' The stacking orientation will be switched (row-stacked becomes column-stacked and vice versa).
#'
#' @param x A multiblock object
#' @return A new multiblock object with transposed matrices
#' @family multiblock functions
#' @export
#' 
#' @examples 
#' mb <- multiblock(list(matrix(1:10, 5, 2), matrix(11:20, 5, 2)))
#' mb_t <- t(mb)
#' is_cstacked(mb_t)
#' is_rstacked(mb_t)
t.multiblock_list <- function(x) {
  out <- lapply(x, t)
  multiblock(out)
}

#' Test if Multiblock Object is Column-Stacked
#'
#' @description
#' Checks if a multiblock object is column-stacked (matrices share row dimension).
#'
#' @param x A multiblock object
#' @return Logical indicating if object is column-stacked
#' @family multiblock functions
#' @seealso [is_rstacked()]
#' @method is_cstacked multiblock_list
#' @export
is_cstacked.multiblock_list <- function(x) {
  inherits(x, "multiblock_c")
}

#' Test if Multiblock Object is Row-Stacked
#'
#' @description
#' Checks if a multiblock object is row-stacked (matrices share column dimension).
#'
#' @param x A multiblock object
#' @return Logical indicating if object is row-stacked
#' @family multiblock functions
#' @seealso [is_cstacked()]
#' @method is_rstacked multiblock_list
#' @export
is_rstacked.multiblock_list <- function(x) {
  inherits(x, "multiblock_r")
}

#' Print Method for Multiblock Objects
#'
#' @description
#' Displays a summary of a multiblock object, including the number of blocks,
#' stacking orientation, and dimensions of each block.
#'
#' @param x A multiblock object
#' @param ... Additional arguments passed to print methods
#' @return Invisibly returns the input object
#' @method print multiblock_list
#' @export
print.multiblock_list <- function(x, ...) {
  cat(crayon::bold(crayon::blue("\n═══ Multiblock Object ═══\n")))
  cat(crayon::bold("\nNumber of blocks: "), crayon::green(length(x)), "\n")
  cat(crayon::bold("Orientation: "), crayon::green(attr(x, "orient")), "\n\n")
  
  if (attr(x, "orient") == "cstacked") {
    cat(crayon::bold("Shared dimension: "), 
        crayon::green(paste(nrow(x[[1]]), "rows")), "\n")
    cat(crayon::bold("Block-specific columns: "), 
        crayon::green(paste(sapply(x, ncol), collapse=", ")), "\n")
    cat(crayon::bold("Total columns: "), 
        crayon::green(sum(sapply(x, ncol))), "\n")
  } else {
    cat(crayon::bold("Shared dimension: "), 
        crayon::green(paste(ncol(x[[1]]), "columns")), "\n")
    cat(crayon::bold("Block-specific rows: "), 
        crayon::green(paste(sapply(x, nrow), collapse=", ")), "\n")
    cat(crayon::bold("Total rows: "), 
        crayon::green(sum(sapply(x, nrow))), "\n")
  }
  
  cat(crayon::bold(blue("\n═══════════════════════\n")))
  invisible(x)
}
